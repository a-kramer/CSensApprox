#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <dirent.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <dlfcn.h>
#include "h5block.h"
#include "ndarray.h"
#include "solution.h"
#include "sensitivity_approximation.h"
#define FALSE 0

#define DEFAULT 1
#define FREE_ON_SUCCESS 1
#define KEEP_ON_SUCCESS 2

typedef struct {
  gsl_vector *u;
  gsl_vector *y0;
  gsl_vector *t;
  char *name;
} simulation_t;

typedef int(*jacp)(double, const double [], double *, void *);


/*This function allocates memory and concatenates two strings in that
  memory. It is used to make function names (this is for loading the model
  functions by name from a shared library). The model function names
  have this pattern: `MODEL_vf`, `MODEL_jac`, `MODEL_jacp`.*/
char* /* string with model_name and suffix (free after loading the function) */
model_function(char *model_name, /* the base name of the model */
 char *suffix) /* suffix, usually `"_vf"` or `"_jac"` */
{
  assert(model_name);
  size_t size=(strlen(model_name)+strlen(suffix)+1);
  assert(size);
  char *f=malloc(sizeof(char)*size);
  assert(f);
  strcat(strcpy(f,model_name),suffix);
  fprintf(stderr,"[%s] «%s»\n",__func__,f); fflush(stderr);
  return f;
}


/* Loads a function from an `.so` file, using `dlsym()`. If dlsym
   fails to find the file `abort()` is called. Optionally, this
   function frees the storage assosiated with the name of the
   function.*/
void *load_or_exit(void *lib, /* file pointer, previously opened via `dlopen()` */
 char *name, /* function to be loaded from file */
 int opt) /* whether to call `free()` on `name` (either: `KEEP_ON_SUCCESS` or `FREE_ON_SUCCESS`). */
{
  assert(lib && name);
  void *symbol=dlsym(lib,name);
  if (symbol) {
    fprintf(stderr,"[%s] «%s» loaded successfully.\n",__func__,name);
    if (opt==FREE_ON_SUCCESS){
      free(name);
    }
  }else{
    fprintf(stderr,"[%s] %s\n",__func__,dlerror());
    abort();
  } 
  return symbol;
}


/* Loads the ODE system from an `.so` file, the file is given by name,
   the returned structure is intended for the `gsl_odeiv2` library of
   solvers. The jacobian dfdx is loaded alongside the right hand side;
   `gsl_odeiv2_system` has no slot for the parameter derivative `dfdp`,
   this is returned as a function pointer instead. */
gsl_odeiv2_system /* the system structure, see gsl documentation. */
load_system(char *model_name, /* the file-name will be constructed from this name, possibly from @link first_so@ */
 size_t n, /* number of state variables */
 double *p, /* default parameter vector */
 jacp *dfdp) /* [output] additional return value: a pointer to the parameter derivative (matrix) function. */
{
  char *so=model_function(model_name,".so");
  char *local_so = model_function("./",so);
  void *lib=dlopen(local_so,RTLD_LAZY);
  free(local_so);
  void *f,*dfdy;
  //dfdp=malloc(sizeof(jacp*));
  char *symbol_name; // symbol name in .so
  if (lib){
    symbol_name=model_function(model_name,"_vf");
    f=load_or_exit(lib,symbol_name,FREE_ON_SUCCESS);
    symbol_name=model_function(model_name,"_jac");
    dfdy=load_or_exit(lib,symbol_name,FREE_ON_SUCCESS);
    symbol_name=model_function(model_name,"_jacp");
    *dfdp=load_or_exit(lib,symbol_name,FREE_ON_SUCCESS);
  } else {
    fprintf(stderr,"[%s] library «%s» could not be loaded: %s\n",__func__,so,dlerror());
    abort();
  }
  gsl_odeiv2_system sys={f,dfdy,n,p};
  free(so);
  fprintf(stderr,"[%s] ode system created.\n",__func__); fflush(stderr);
  return sys;
}


/* A command line option has the pattern `"-o"` or `"--option"`. This
   function checks whether a given string `value` is the name of a
   short option `s` or long option `l`.*/
int /* returns `1` if the given string is the specified option, `0` otherwise */
option_is(char *s, /* short option name, e.g.: "-t", can be NULL */
 char *l, /* long option name, e.g.: "--time", can also be NULL */
 char *value) /* a string for comaprison with s and l*/
{
  int match=FALSE;
  assert(value);
  if (s) match=(strcmp(value,s)==0);
  if (l) match|=(strcmp(value,l)==0);
  return match;
}


/* loads the simulation instructions from an hdf5 file, contained in
   an hdf5 group (g_id). Simulation instructions are: the initial
   value, the input parameters, and time scale of the problem. This
   function has the interface the is required for `H5Literate()`. */
herr_t /* a nonnegative value is interpreted as success by the hdf5 library. */
load_one_sim(hid_t g_id, /* the group identifier, this group holds the simulation instructions*/
 const char *name, /* the name of the group */
 const H5L_info_t *info, /* hdf5 specific object (not used) */
 void *op_data) /* [OUT] an array of simulation structures. */
{
  gsl_vector *u = h5_to_gsl(g_id,name,"input");
  gsl_vector *y0 = h5_to_gsl(g_id,name,"InitialValue");
  gsl_vector *t = h5_to_gsl(g_id,name,"time");
  gsl_vector_int *index = h5_to_gsl_int(g_id,name,"index");
  simulation_t *sim = op_data;
  int i=gsl_vector_int_get(index,0);
  fprintf(stderr,"[%s] processing «%s» (index %i)\n",__func__,name,i);
  assert(u && y0 && t);
  sim[i].u=u;
  sim[i].y0=y0;
  sim[i].t=t;
  sim[i].name=strdup(name);
  return 0;
}


/* Loads all simulation instructions (experiment
   descriptions/protocols) from an hdf5 file using `H5Literate()` */
simulation_t* /* returns `N` simulations as a structure array (allocated via `malloc()`) */
sim_from_h5(hid_t g_id, /* group id, location id of all simulation DATASETS */
 hsize_t *N) /* [OUT] additionally returns the number of specified simulations `N`. */
{
  herr_t h5ec=H5Gget_num_objs(g_id,N);
  assert(h5ec>=0);
  fprintf(stderr,"[%s] number of data sets: %lli\n",__func__,N[0]); fflush(stderr);
  simulation_t *sim=malloc(sizeof(simulation_t)*N[0]);
  hsize_t idx=0;
  h5ec=H5Literate(g_id, H5_INDEX_NAME, H5_ITER_NATIVE, &idx, load_one_sim, sim);
  assert(h5ec>=0);
  return sim;
}


/* If no simulation time span was specified, it is inferred from the
   measurement times `t` in the simulation instructions. The `tspan`
   vector has the same three components as GNU Octave uses to specify
   a range of numbers: `initial` `increment` `final` (e.g.: `"0:0.1:7"`).
   If tspan is `[0,0,0]`, then it is corrected to 100 points between `0` and
   `1.3*max(t)` */
void fix_tspan_if_necessary(gsl_vector *t, /* measurement time vector (from the data file) */
 ndarray *tspan) /* 3-element tspan vector to fix */
{
  assert(tspan);
  int i=0;
  double tf=gsl_vector_max(t)*1.3; // a reasonable stop time [default value]
  double *a,*b=NULL,*c;

  a=ndarray_ptr(tspan,&i);
  if (tspan->size[0]>3) {
    i++;
    b=ndarray_ptr(tspan,&i);
  }
  i++;
  c=ndarray_ptr(tspan,&i);
  
  if (*a==*c) *c=tf;
  assert(*c>*a);
  if (b && *b==0.0) *b=1e-2 * (*c - *a);
  fprintf(stderr,"[%s] simulating in t: [%g,%g]",__func__,*a,*c);
  if (b) fprintf(stderr," with an increment of %g\n",*b);
}



void save_solution(gsl_odeiv2_system sys, jacp dfdp, solution_t *solution, gsl_vector *y, gsl_vector* par, int *index, double t){
  double *Jy,*Jp,*y_ptr,*f;
  int d=y->size;
  double dfdt[d];
  y_ptr = ndarray_ptr(solution->y,&(index[1]));
  f = ndarray_ptr(solution->f,&(index[1]));
  Jy = ndarray_ptr(solution->Jy,index);
  Jp = ndarray_ptr(solution->Jp,index);
  sys.function(t, y->data, f, par->data);
  sys.jacobian(t, y->data, Jy, dfdt, par->data);
  dfdp(t, y->data, Jp, par->data);
  
  memcpy(y_ptr,y->data,sizeof(double)*d);
  *ndarray_ptr(solution->t,&(index[2]))=t;
}

void stop_if_unsuccessful(int status, char *name, double t, double tf){
  switch (status){
  case GSL_SUCCESS:
    return;
  case GSL_EMAXITER:
    fprintf(stderr,"[%s] simulation %s, time %g: maximum number of steps reached.\n",__func__,name,t);
    fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
    break;
  case GSL_ENOPROG:
    fprintf(stderr,"[%s] simulation %s, time %g: step size dropeed below set minimum.\n",__func__,name,t);
    fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
    break;
  case GSL_EBADFUNC:
    fprintf(stderr,"[%s] simulation %s: bad function.\n",__func__,name);
    fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
    break;
  default:
    fprintf(stderr,"[%s] unhandled error code: %#x\n",__func__,status);
    abort();
  }      
  printf("\n\n");
}

void write_solution(hid_t h5f, char *name, solution_t* solution){
  hid_t g_id;
  if (h5f) {
    g_id = H5Gcreate2(h5f, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ndarray_to_h5(solution->y,g_id,"state");
    ndarray_to_h5(solution->f,g_id,"f");
    ndarray_to_h5(solution->t,g_id,"time");
    ndarray_to_h5(solution->status,g_id,"status");
    ndarray_to_h5(solution->Jy,g_id,"jac");
    ndarray_to_h5(solution->Jp,g_id,"jacp");
    ndarray_to_h5(solution->Sy,g_id,"sensitivity");
    ndarray_to_h5(solution->PHIf,g_id,"transition_matrix_forward");
    ndarray_to_h5(solution->PHIb,g_id,"transition_matrix_backward");
    
    H5Gclose(g_id);
  }
}



/* Intergrates the system `sys` using the specified `driver` and
   simulation instructions `sim` (an array of structs, one element per
   simulation). The results are saved to an hdf5 file and also printed
   to standard output. */
solution_t** /* a structure that contains the trajectories y(t) as well as jacobians Jy, Jp and the y-sensitivity (dy/dp).*/
simulate(gsl_odeiv2_system sys, /* the system to integrate */
 gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
 simulation_t *sim, /*`N` simulation instructions (e.g.: initial value y0)*/
 hsize_t N, /* number of simulations to perform, length of `sim` */
 ndarray *tspan, /* time span vector: initial time, increment, final time */
 gsl_vector *u, /* input vector (a pointer that can be used to change `sys`)*/ 
 gsl_vector *par, /* parameters of the model, a pointer that can be used to change `sys`.*/
 jacp dfdp, /* a function that returns the parameter derivative of the model's right hand side function.*/
 hid_t h5f) /* an hdf5 file opened for writing.*/
{
  assert(driver && sim && N>0);
  assert(tspan);
  assert(sim);
  size_t d=sim[0].y0->size;
  size_t p=par->size;
  gsl_vector *y=gsl_vector_alloc(d);
  size_t i,j,k;
  solution_t **solution=malloc(sizeof(solution_t*)*N);
  int index[3]={0,0,0};
  double t,tf,dt;
  size_t n;
  int status;
  for (i=0;i<N;i++){
    assert(sim[i].name && strlen(sim[i].name));
    gsl_odeiv2_driver_reset(driver);
    // fix possible issues with simulation time specification:
    fix_tspan_if_necessary(sim[i].t,tspan);    
    gsl_vector_memcpy(y,sim[i].y0);
    //gsl_vector_add_constant(y,1e-5);
    gsl_vector_memcpy(u,sim[i].u);
     t=tspan->value[0];
    dt=tspan->value[1];
    // print table header:
    printf("# Simulation %li: \n",i);
    printf("#t\t");
    for (k=0;k<y->size;k++)
      printf("y%li\t",k);
    printf("\n");
    n=(size_t) ((tspan->value[2]-tspan->value[0])/tspan->value[1]);
    solution[i]=solution_alloc(d,p,n);
    for (j=0;j<n;j++){
      tf=t+dt;
      status=gsl_odeiv2_driver_apply(driver, &t, tf, y->data);
      //report any error codes to the user
      stop_if_unsuccessful(status,sim[i].name,t,tf);
      printf("%g\t",t);
      for (k=0;k<y->size;k++)
	printf("%g\t",gsl_vector_get(y,k));
      printf("\n");
      index[2]=j;
      save_solution(sys, dfdp, solution[i], y, par, index, t);
    }
    sensitivity_approximation(solution[i]);
    write_solution(h5f,sim[i].name,solution[i]);
  }
  return solution;
}

/* Intergrates the system `sys` using the specified `driver`, but
   using the solvers steps, ignoring tspan[1]. The simulation
   instructions are stored in `sim` (an array of structs, one element
   per simulation). The results are saved to an hdf5 file and also
   printed to standard output. */
solution_t** /* a structure that contains the trajectories y(t) as well as jacobians Jy, Jp and the y-sensitivity (dy/dp).*/
simulate_evolve(gsl_odeiv2_system sys, /* the system to integrate */
 gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
 simulation_t *sim, /*`N` simulation instructions (e.g.: initial value y0)*/
 hsize_t N, /* number of simulations to perform, length of `sim` */
 ndarray *tspan, /* time span vector: initial time, final time */
 gsl_vector *u, /* input vector (a pointer that can be used to change `sys`)*/ 
 gsl_vector *par, /* parameters of the model, a pointer that can be used to change `sys`.*/
 jacp dfdp, /* a function that returns the parameter derivative of the model's right hand side function.*/
 hid_t h5f) /* an hdf5 file opened for writing.*/
{
  assert(driver && sim && N>0);
  assert(tspan);
  assert(sim);
  size_t d=sim[0].y0->size;
  size_t p=par->size;
  gsl_vector *y=gsl_vector_alloc(d);
  double h=1e-2;
  size_t i,k;
  solution_t **solution=malloc(sizeof(solution_t*)*N);
  int index[3]={0,0,0};
  double t,tf;
  size_t n,n_max=100;
  int status;
  for (i=0;i<N;i++){// simulations
    assert(sim[i].name && strlen(sim[i].name));
    gsl_odeiv2_driver_reset(driver);
    // fix possible issues with simulation time specification:
    fix_tspan_if_necessary(sim[i].t,tspan);
    gsl_vector_memcpy(y,sim[i].y0);
    //gsl_vector_add_constant(y,1e-5);
    gsl_vector_memcpy(u,sim[i].u);
    t=tspan->value[0];
    tf=tspan->value[1];
    // print table header:
    printf("# Simulation %li: \n",i);
    printf("#t\t");
    for (k=0;k<y->size;k++) printf("y%li\t",k);
    printf("\n");
    // estimate number of points:
    n=0;
    solution[i]=solution_alloc(d,p,n_max);
    while (t<tf){
      if (n==n_max){
	n_max+=100;
	solution_resize(solution[i],n_max);
      }      
      printf("%g\t",t);
      for (k=0;k<y->size;k++) printf("%g\t",gsl_vector_get(y,k));
      printf("\n");
      index[2]=n;
      save_solution(sys, dfdp, solution[i], y, par, index, t);
      status = gsl_odeiv2_evolve_apply
	(driver->e, /* evolve object        */
	 driver->c, /* control object       */
	 driver->s, /* stepper object       */
	 &sys,      /* system               */
	 &t,        /* time to update       */
	 tf,        /* direction            */
	 &h,        /* step size suggestion */ 
	 y->data);
      stop_if_unsuccessful(status,sim[i].name,t,tf);
      n++;
    }
    solution_resize(solution[i],n);
    sensitivity_approximation(solution[i]);
    write_solution(h5f,sim[i].name,solution[i]);
  }
  return solution;
}

/*Finds the first file ending in `.so` in the current working
  directory. The suffix `.so` is stripped from the file name, the
  prefix is assumed to be a the model's name.  */
char* /* the name of the found model, NULL if no file was found.*/
first_so(){
  char *name = NULL;
  char *p;
  DIR *dp = opendir (".");
  assert(dp);
  struct dirent *d;
  while ((d=readdir (dp))!=NULL){
    p=strrchr(d->d_name,'.');
    if (d->d_type==DT_REG && p && strcmp(p,".so")==0){
      p[0]='\0';
      name=strdup(d->d_name);
      p[0]='.';
      break;
    }
  }
  closedir (dp);
  if (name){
    fprintf(stderr,"[%s] automatically determined model name: «%s».\n",__func__,name);
  } else {
    fprintf(stderr,"[%s] no model name specified: please use the --model option.\n",__func__);
    abort();
  }
  return name;
}

/* Interprets a string as a range specification from three values:
   "initial increment final", the values can be separated by spaces or
   colons `:`. */
ndarray* /* a three element array of `double`s */
read_tspan(char *val)/*a string of the form "a:b:c" or "a b c". */{
  int n=0,j;
  char *s,*p;
  int size=3;
  ndarray *t=ndarray_alloc(1,&size);
  if (val){
    s=val;
    p=val;
    for (j=0;j<3;j++){
      if (p[0]==':') p++;
      s=p;
      *ndarray_ptr(t,&j)=strtod(s,&p);
      if (s==p) break;
      else n++;
    }
    ndarray_resize(t,n);
  }
  return t;
}

void test_evaluation(gsl_odeiv2_system sys, jacp dfdp, gsl_vector *y0, gsl_vector *par){
  //test evaluation:
  size_t d=sys.dimension;
  int jacp_size=d*par->size;
  double *Jy=malloc(sizeof(double)*(d*d));
  double *Jp=malloc(sizeof(double)*jacp_size);
  double *f=malloc(sizeof(double)*d);
  int i;
  // f
  sys.function(0,y0->data,f,par->data);
  fprintf(stderr,"[%s] test evaluation of right hand side function (f):\n",__func__);
  for (i=0;i<d;i++) fprintf(stderr,"%g ",f[i]);
  fprintf(stderr,"\n\n");
  fflush(stderr);

  //jacp
  dfdp(0,y0->data,Jp,par->data);
  fprintf(stderr,"[%s] test evaluation of jacp (df/dp):\n",__func__);
  for (i=0;i<jacp_size;i++) fprintf(stderr,"%g ",Jp[i]);
  fprintf(stderr,"\n (that was a flat %li × %li matrix)\n",d,par->size);
  fflush(stderr);
  // jac

  sys.jacobian(0,y0->data,Jy,f,par->data);
  fprintf(stderr,"[%s] test evaluation of jacobian (df/dy):\n",__func__);
  for (i=0;i<d*d;i++) fprintf(stderr,"%g ",Jy[i]);
  fprintf(stderr,"\n (that was a flat %li × %li matrix)\n",d,d);
  fflush(stderr);
  free(Jy);
  free(Jp);
  free(f);


}

/* This prgram loads an ODE model, specified for the `gsl_odeiv2`
   library. The initial value problems are specified in an hdf5 file,
   intended for use in systems biology applications. So, some of the
   terminology in the expected hdf5 _data_ file is vaguely related to
   biological systems. All command line arguments are optional and
   names of files are guessed, based on the contents of the current
   working directoy. The hdf5 fiel is expected to have a group called
   "data", this group shall contain hdf5 DATASETS with ATTRIBUTES that
   describe initial value problems suitable to replicate these
   datasets: InitialValue, time, and index. Currenty, the model is
   parameterized using the value of the DATASET *mu*, in the GROUP
   called *prior*. 

   The possible command line options are documented in the [README.md](../README.md) .
*/
int /* `EXIT_SUCESS` if all files are found and integration succeeds, default `abort()` signal otherwise.*/
main(int argc, char *argv[]){
  int i=0;
  char *model_name=NULL, *h5file=NULL;
  double abs_tol=1e-6,rel_tol=1e-5,h=1e-3;
  ndarray *t=NULL;
  ndarray_test();
  for (i=1;i<argc;i++){
    fprintf(stderr,"[%s] %s=%s\n",__func__,argv[i],argv[i+1]); fflush(stderr);
    if (option_is("-m","--model",argv[i])){
      i++;
      model_name=strdup(argv[i]); //malloc(sizeof(char)*(strlen(argv[i])+1));
      //strcpy(model_name,argv[i]);
    } else if (option_is("-d","--data",argv[i])){
      i++;
      h5file=strdup(argv[i]); //malloc(sizeof(char)*(strlen(argv[i])+1));
      //strcpy(h5file,argv[i]);
    } else if (option_is("-a","--abs-tol",argv[i])){
      i++;
      abs_tol=strtod(argv[i],NULL);
    } else if (option_is("-r","--rel-tol",argv[i])){
      i++;
      rel_tol=strtod(argv[i],NULL);
    } else if (option_is("-h","--step-size",argv[i])){
      i++;
      h=strtod(argv[i],NULL);
    } else if (option_is("-t","--sim-time",argv[i])){
      i++;
      t=read_tspan(argv[i]);
    }else {
      fprintf(stderr,"[%s] unknown option «%s»\n",__func__,argv[i]);
      exit(1);
    }
  }
  if (!t) t=read_tspan("0 1");
  if (!model_name) model_name=first_so();
  assert(model_name);
  
  if (!h5file) h5file=model_function(model_name,".h5");
  assert(h5file);
  hid_t h5f_id=H5Fopen(h5file,H5F_ACC_RDONLY,H5P_DEFAULT);
  assert(h5f_id);
  hid_t prior=H5Gopen2(h5f_id,"prior",H5P_DEFAULT);
  gsl_vector *mu=h5_to_gsl(prior,"mu",NULL);
 
  fprintf(stderr,"[%s] prior: ",__func__);
  for (i=0;i<mu->size;i++) fprintf(stderr,"%g ",gsl_vector_get(mu,i));
  fprintf(stderr,"\n");
  
  hsize_t N;
  hid_t data=H5Gopen2(h5f_id,"data",H5P_DEFAULT);
  fprintf(stderr,"[%s] data id: %ld\n",__func__,data);
  simulation_t *sim = sim_from_h5(data,&N);
  size_t d=sim[0].y0->size; // get from initial conditions

  size_t nu=sim[0].u->size;
  gsl_vector *par=gsl_vector_alloc(mu->size + nu);
  gsl_vector_view p=gsl_vector_subvector(par,0,mu->size);
  gsl_vector_view u=gsl_vector_subvector(par,mu->size,nu);
  // initialize:
  gsl_vector_memcpy(&(p.vector),mu);
  gsl_vector_memcpy(&(u.vector),sim[0].u);
  
  for (i=0;i<mu->size;i++) par->data[i]=exp(par->data[i]);

  // load system from file and test it
  jacp dfdp;
  gsl_odeiv2_system sys = load_system(model_name, d, par->data, &dfdp);
  fprintf
    (stderr,
     "[%s] ode system dim: %li and %li parameters (%li inputs) with parameters:\n",
     __func__,
     sys.dimension,
     par->size,
     nu);
  for (i=0;i<mu->size;i++) fprintf(stderr,"%g ",((double*) sys.params)[i]);
  fprintf(stderr,"(that is exp(mu))\n");
  
  test_evaluation(sys,dfdp,sim[0].y0,par);

  const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;
  gsl_odeiv2_driver *driver=gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
  fprintf(stderr,"[%s] driver allocated with tolearances: abs_tol=%g and rel_tol=%g\n",
	 __func__,abs_tol,rel_tol);
  char *h5out_name = model_function(model_name,"_out.h5");
  fprintf(stderr,"[%s] output file: %s\n",__func__,h5out_name); 
  hid_t h5f = H5Fcreate(h5out_name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  assert(h5f);

  // the most CPU work happens here:
  if (t->size[0]==2){
    solution_t **solution=simulate_evolve(sys,driver,sim,N,t,&(u.vector),par,dfdp,h5f);
  } else if (t->size[0]==3){
    solution_t **solution=simulate(sys,driver,sim,N,t,&(u.vector),par,dfdp,h5f);
  }
  gsl_odeiv2_driver_free(driver);
  H5Fclose(h5f);
  return EXIT_SUCCESS;
}
