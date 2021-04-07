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
#include <math.h>
#include <dlfcn.h>
#include "h5block.h"
#include "ndarray.h"
#include "solution.h"
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


char* model_function(char *model_name, char *suffix){
  assert(model_name);
  size_t size=(strlen(model_name)+strlen(suffix)+1);
  assert(size);
  char *f=malloc(sizeof(char)*size);
  assert(f);
  //*((char *) mempcpy(mempcpy(f,model_name,strlen(model_name)),suffix,strlen(suffix)))='\0';
  strcat(strcpy(f,model_name),suffix);
  fprintf(stderr,"[%s] «%s»\n",__func__,f); fflush(stderr);
  return f;
}

void *load_or_exit(void *lib, char *name, int opt){
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

gsl_odeiv2_system load_system(char *model_name, size_t n, double *p, jacp *dfdp){
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

int option_is(char *s, char *l, char *value){
  int match=FALSE;
  assert(value);
  if (s) match=(strcmp(value,s)==0);
  if (l) match|=(strcmp(value,l)==0);
  return match;
}

herr_t load_one_sim(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data){
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

simulation_t* sim_from_h5(hid_t g_id, hsize_t *N){
  herr_t h5ec=H5Gget_num_objs(g_id,N);
  assert(h5ec>=0);
  fprintf(stderr,"[%s] number of data sets: %lli\n",__func__,N[0]); fflush(stderr);
  simulation_t *sim=malloc(sizeof(simulation_t)*N[0]);
  hsize_t idx=0;
  h5ec=H5Literate(g_id, H5_INDEX_NAME, H5_ITER_NATIVE, &idx, load_one_sim, sim);
  assert(h5ec>=0);
  return sim;
}

void fix_tspan_if_necessary(gsl_vector *t, double *tspan){
  assert(tspan);
  double tf=gsl_vector_max(t)*1.3; // a reasonable stop time [default value]
  if (tspan[2]==tspan[0]) tspan[2]=tf;
  assert(tspan[2]>tspan[0]);
  if (tspan[1]==0.0) tspan[1]=1e-2*(tspan[2]-tspan[0]);
  fprintf(stderr,"[%s] simulating in t: [%g,%g] with an increment of %g\n",__func__,tspan[0],tspan[2],tspan[1]);

}


gsl_matrix* transition_matrix(gsl_matrix *Ji, gsl_matrix *Jf, double ti, double tf){
  double s=0.5*(tf-ti);
  double b=0.0;
  size_t ny=Ji->size1;
  gsl_matrix *identity=gsl_matrix_alloc(ny,ny);
  gsl_matrix *PHI=gsl_matrix_alloc(ny,ny);
  gsl_matrix_set_identity(PHI);
  gsl_matrix_set_identity(identity);
  int i,n=3;
  for (i=0;i<n;i++){
    // PHI = s*Jf*PHI + b*PHI + I:
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, s, Jf, PHI, b, PHI);
    gsl_matrix_add(PHI,identity);
  }
  gsl_matrix_free(identity);
  return PHI;
}
     
void sensitivity_approximation(solution_t *solution){
  assert(solution && solution->Jp && solution->Jy && solution->t);

  int np = solution->Jp->size[0];
  int ny = solution->Jp->size[1];
  int nt = solution->t->size;
  int index[3]={0,0,0};
  int j;
  double tf, ti;
  gsl_matrix_view Sf,Si,Ji,Jf,B;
  gsl_matrix *PHI_fwd, *PHI_bwd;
  gsl_matrix *S_temp=gsl_matrix_alloc(ny,np);
  assert(S_temp);

  for (j=1;j<nt;j++){
    index[2]=j; // time index
    Sf = gsl_matrix_view_array(ndarray_ptr(solution->Sy,index),ny,np);
    Jf=gsl_matrix_view_array(ndarray_ptr(solution->Jy,index),ny,ny);
    B=gsl_matrix_view_array(ndarray_ptr(solution->Jp,index),ny,np);
    tf=gsl_vector_get(solution->t,j);
    
    index[2]=j-1; // previous time index
    Si = gsl_matrix_view_array(ndarray_ptr(solution->Sy,index),ny,np);
    Ji=gsl_matrix_view_array(ndarray_ptr(solution->Jy,index),ny,ny);
    ti=gsl_vector_get(solution->t,j-1);

    PHI_fwd = transition_matrix(&Ji.matrix,&Jf.matrix,ti,tf);
    PHI_bwd = transition_matrix(&Ji.matrix,&Jf.matrix,tf,ti);
    // Sf = PHI(k+1,k) * (0.5*(PHI(k,k+1)*B(k) + B(k)) + Si )
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, PHI_bwd, &B.matrix, 0.0, S_temp);
    gsl_matrix_add(S_temp,&B.matrix);
    gsl_matrix_scale(S_temp,0.5);
    gsl_matrix_add(S_temp,&Si.matrix);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, PHI_fwd, S_temp, 0.0, &Sf.matrix);
  }
  gsl_matrix_free(S_temp);
}


solution_t** simulate(gsl_odeiv2_system sys, gsl_odeiv2_driver* driver, simulation_t *sim, hsize_t N, double *tspan, gsl_vector *u, gsl_vector *par, jacp dfdp, hid_t h5f){
  
  assert(driver && sim && N>0);
  assert(tspan);
  assert(sim);
  size_t d=sim[0].y0->size;
  size_t p=par->size;
  gsl_vector *y=gsl_vector_alloc(d);
  gsl_vector_view y_row;
  size_t i,j,k;
  solution_t **solution=malloc(sizeof(solution_t*)*N);
  int index[3]={0,0,0};
  double t,tf,dt;
  double *Jy,*Jp;
  double *dfdt=malloc(sizeof(double)*d);
  size_t n;
  hid_t g_id;
  int status;
  for (i=0;i<N;i++){
    assert(sim[i].name && strlen(sim[i].name));
    if (h5f) g_id = H5Gcreate2(h5f, sim[i].name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    gsl_odeiv2_driver_reset(driver);
    // fix possible issues with simulation time specification:
    fix_tspan_if_necessary(sim[i].t,tspan);    
    gsl_vector_memcpy(y,sim[i].y0);
    //gsl_vector_add_constant(y,1e-5);
    gsl_vector_memcpy(u,sim[i].u);
     t=tspan[0];
    dt=tspan[1];
    // print table header:
    printf("# Simulation %li: \n",i);
    printf("#t\t");
    for (k=0;k<y->size;k++)
      printf("y%li\t",k);
    printf("\n");
    n=(size_t) ((tspan[2]-tspan[0])/tspan[1]);
    solution[i]=solution_alloc(d,p,n);
    for (j=0;j<n;j++){
      tf=t+dt;
      status=gsl_odeiv2_driver_apply(driver, &t, tf, y->data);
      //report any error codes to the user
      switch (status){
      case GSL_EMAXITER:
	fprintf(stderr,"[%s] simulation %li, time_point %li: maximum number of steps reached.\n",__func__,i,j);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_ENOPROG:
	fprintf(stderr,"[%s] simulation %li, time_point %li: step size dropeed below set minimum.\n",__func__,i,j);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_EBADFUNC:
	fprintf(stderr,"[%s] simulation %li, time_point %li: bad function.\n",__func__,i,j);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_SUCCESS:
	printf("%g\t",t);
	for (k=0;k<y->size;k++)
	  printf("%g\t",gsl_vector_get(y,k));
	printf("\n");
	index[2]=j;
	Jy = ndarray_ptr(solution[i]->Jy,index);
	Jp = ndarray_ptr(solution[i]->Jp,index);
	sys.jacobian(t, y->data, Jy, dfdt, par->data);
	dfdp(t,y->data,Jp,par->data);
	y_row = gsl_matrix_row(solution[i]->y,j);
	gsl_vector_memcpy(&y_row.vector,y);
	gsl_vector_set(solution[i]->t,j,t);
	break;
      default:
	fprintf(stderr,"[%s] unhandled error code: %#x\n",__func__,status);
	abort();
      }
      if (status!=GSL_SUCCESS) break;
    }
    sensitivity_approximation(solution[i]);
    if (h5f){
      gsl_matrix_to_h5(solution[i]->y,g_id,"state");
      gsl_vector_to_h5(solution[i]->t,g_id,"time",NULL);
      ndarray_to_h5(solution[i]->Jy,g_id,"jac");
      ndarray_to_h5(solution[i]->Jp,g_id,"jacp");
      ndarray_to_h5(solution[i]->Sy,g_id,"sensitivity");
      H5Gclose(g_id);
    }    
    printf("\n\n");
  }
  return solution;
}

char *first_so(){
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

double* read_tspan(char *val){
  int j;
  char *s,*p;
  double *t=calloc(3,sizeof(double));
  if (val){
    s=val;
    p=val;
    for (j=0;j<3;j++){
      if (p[0]==':') p++;
      s=p;
      t[j]=strtod(s,&p);
      if (s==p) break;
    }
  }
  return t;
}

int main(int argc, char *argv[]){
  int i=0;
  char *model_name=NULL, *h5file=NULL;
  double abs_tol=1e-6,rel_tol=1e-5,h=1e-3;
  double *t=NULL;
  ndarray_test();
  for (i=1;i<argc;i++){
    fprintf(stderr,"[%s] %s=%s\n",__func__,argv[i],argv[i+1]);
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
  if (!t) t=read_tspan("0 0 0");
  fprintf(stderr,"[%s] t: %g:%g:%g\n",__func__,t[0],t[1],t[2]);

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
  gsl_vector_memcpy(&(p.vector),mu);
  for (i=0;i<mu->size;i++) par->data[i]=exp(par->data[i]);
  jacp dfdp;
  gsl_odeiv2_system sys = load_system(model_name, d, par->data, &dfdp);
  //test evaluation:
  int jacp_size=d*par->size;
  double *J=malloc(sizeof(double)*jacp_size);
  gsl_vector_memcpy(&(u.vector),sim[0].u);
  dfdp(0,sim[0].y0->data,J,par->data);
  fprintf(stderr,"[%s] test evaluation of jacp\n",__func__);
  for (i=0;i<jacp_size;i++) fprintf(stderr,"%g ",J[i]);
  fprintf(stderr,"\n\n");
  const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;
  gsl_odeiv2_driver *driver=gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
  char *h5out_name = model_function(model_name,"_out.h5");
  hid_t h5f = H5Fcreate(h5out_name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  assert(h5f);
  // the most CPU work happens here:

  solution_t **solution=simulate(sys,driver,sim,N,t,&(u.vector),par,dfdp,h5f);
  gsl_odeiv2_driver_free(driver);
  H5Fclose(h5f);
  return EXIT_SUCCESS;
}
