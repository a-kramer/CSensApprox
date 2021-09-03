

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



solution_t **solution_from_gsl(){
    jacp dfdp;
    gsl_odeiv2_system sys = load_system(model_name, d, par->data, &dfdp);
    fprintf
      (stderrn,
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
}
