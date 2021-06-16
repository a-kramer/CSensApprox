#include "sensitivity_approximation.h"

/* this function calculates two projections, one projection to state-space that is in steady state; the second    projection to non-steady state variables: 
   let y'=f(t,y), and
   y=[v;w]; where v is in steady state and w is not. More generally
   v=Py
   w=Qy
   where P and Q are projection matrices (so, v and w don't have to be in that order)
   then P*f(t,y) == 0
   and  Q*f(t,y) != 0
   but instead of matrices, we just keep track of the index-structure of the projection:
   v[i]=y[p[i]], where p[i] is a list of indices.
    */
int* not_in_steady_state(solution_t *solution, int** pptr, double abs_tol, double rel_tol){
  ndarray *y=solution->y;
  ndarray *f=solution->f;
  int rank=y->rank;
  assert(rank==2);
  int index[2]={0,0};
  int i,j,k;
  int nt=y->size[1];
  int ny=y->size[0];
  int *s=malloc(sizeof(int)*nt);
  int *p=malloc(sizeof(int)*(nt*ny));
  //ndarray_print(f,"f");
  //ndarray_print(y,"y");
  for (j=0;j<nt;j++){
    s[j]=0;
    index[1]=j;
    k=ny-1;
    for (i=0;i<ny;i++){
      index[0]=i;
      if (fabs(ndarray_value(f,index)) > abs_tol + fabs(ndarray_value(y,index))*rel_tol){
	p[ny*j+s[j]]=i;
	s[j]++;
      } else {
	p[ny*j + k]=i;
	k--;
      }
    }
  }
  *pptr=p;
  return s;
}

/* a state projection to the index set p, of length n*/
state_t* /**/
state_projection
(int n,
 int *p,
 const double t,
 const gsl_matrix *A,
 const gsl_matrix *B,
 const gsl_matrix *S)
{
  assert(A);
  assert(B);
  assert(S);
  state_t *s=NULL;
  int np=B->size2;
  int i,j;
  double d;
  if (n){
    s=malloc(sizeof(state_t));
    s->t=t;
    s->dfdy=gsl_matrix_alloc(n,n);
    s->dfdp=gsl_matrix_alloc(n,np);
    s->dydp=gsl_matrix_alloc(n,np);
    for (i=0;i<n;i++){
      for (j=0;j<n;j++){
	d=gsl_matrix_get(A,p[i],p[j]);
	gsl_matrix_set(s->dfdy,i,j,d);
      }
      for (j=0;j<np;j++){
	d=gsl_matrix_get(B,p[i],j);
	gsl_matrix_set(s->dfdp,i,j,d);
      }
      for (j=0;j<np;j++){
	d=gsl_matrix_get(S,p[i],j);
	gsl_matrix_set(s->dydp,i,j,d);
      }
    }
  }
  return s;
}

void state_projection_free(state_t *s){
    free(s->dfdy);
    free(s->dfdp);
    free(s->dydp);
    free(s);
}

/* Calculates the transition matrix using a truncated Peano Baker
   series, and 1-step trapezoidal approximation of integrals: `tf-ti`
   needs to be small.*/
void 
transition_matrix(gsl_matrix *Ji, /* the jacobian at t=ti */
 gsl_matrix *Jf, /* the jacobian at t=tf */
 double ti, /* initial time of the interval (left bondary) */
 double tf, /* final time of the interval (right boundary) */
 gsl_matrix *phi) /* [OUT] PHI(tf,ti): the transition matrix between ti and tf*/
{
  double s=0.5*(tf-ti);
  size_t ny=Ji->size1;
  gsl_matrix *I_k=gsl_matrix_alloc(ny,ny);   // I_{ k }(tf,ti)
  gsl_matrix *I_k_plus_1=gsl_matrix_alloc(ny,ny);  // I_{k+1}(tf,ti)

  gsl_matrix_set_identity(phi); // I_0
  int i,n=2;
  // I_1 is:
  gsl_matrix_memcpy(I_k,Jf);
  gsl_matrix_add(I_k,Ji);
  gsl_matrix_scale(I_k,s);
  
  for (i=0;i<n;i++){
    gsl_matrix_add(phi,I_k);
    // C = s*A*B + b*C   [dgemm] s is DeltaT*0.5 and b is 0
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, s, Jf, I_k, 0.0, I_k_plus_1);
    gsl_matrix_memcpy(I_k,I_k_plus_1); 
  }
  gsl_matrix_free(I_k);
  gsl_matrix_free(I_k_plus_1);
}

/* Calculates the transition matrix using a truncated Peano Baker
   series, and 1-step trapezoidal approximation of integrals: `tf-ti`
   needs to be small.*/
void transition_matrix_v2(gsl_matrix *Ji, /* the jacobian at t=ti */
 gsl_matrix *Jf, /* the jacobian at t=tf */
 double ti, /* initial time of the interval (left bondary) */
 double tf,/* final time of the interval (right boundary) */
 gsl_matrix *phi) /* return buffer */
{
  double s=0.5*(tf-ti);
  size_t ny=Ji->size1;
  int i,n=3; 
  gsl_matrix *V=gsl_matrix_alloc(ny,ny);
  gsl_matrix *W=gsl_matrix_alloc(ny,ny);
  gsl_matrix *I0=gsl_matrix_alloc(ny,ny);  // I0(tf;ti) = identity;
  gsl_matrix *I1=gsl_matrix_alloc(ny,ny);  // I1(tf;ti) = 0.5*(tf-ti)*(Jf+Ji)
  gsl_matrix_set_identity(phi);
  gsl_matrix_set_identity(I0);
  // I1
  gsl_matrix_memcpy(I1,Jf);
  gsl_matrix_add(I1,Ji);
  gsl_matrix_scale(I1,s);
  
  gsl_matrix_set_identity(W);
  for (i=0;i<n;i++){
    gsl_matrix_set_identity(V);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, s, W, Jf, 1.0, V);
    gsl_matrix_memcpy(W,V);
  }
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, W, I1, 1.0, phi);

  gsl_matrix_free(I0);
  gsl_matrix_free(I1);
  gsl_matrix_free(V);
  gsl_matrix_free(W);
}

expm_work_t* work_mem_alloc(int ny, int np){
  expm_work_t *w=malloc(sizeof(expm_work_t));
  w->A=gsl_matrix_alloc(ny,ny);
  w->eA=gsl_matrix_alloc(ny,ny);
  w->LU=gsl_matrix_alloc(ny,ny);
  w->prm=gsl_permutation_alloc(ny);
  w->A_B=gsl_matrix_alloc(ny,np);
  w->S=gsl_matrix_alloc(ny,np);
  return w;
}

void work_mem_free(expm_work_t *w){
  gsl_matrix_free(w->eA);
  gsl_matrix_free(w->A);
  gsl_matrix_free(w->LU);
  gsl_matrix_free(w->A_B);
  gsl_permutation_free(w->prm);
  free(w);
}

void
near_steady_state_approximation
(gsl_matrix *Sf,
 gsl_matrix *Si,
 gsl_matrix *Ji,
 gsl_matrix *Bi,
 double tf,
 double ti,
 expm_work_t *w)
{
  int i;
  //int ny=Bi->size1;
  int np=Bi->size2;
  gsl_vector_view x,b,diag;
  int status;

  gsl_matrix_memcpy(w->A,Ji);
  diag=gsl_matrix_diagonal(w->A);
  gsl_vector_add_constant(&diag.vector,1e-10);
  gsl_matrix_memcpy(w->LU,w->A);
  status=gsl_linalg_LU_decomp(w->LU, w->prm, &(w->sign));
  assert(status==GSL_SUCCESS);
  
  gsl_matrix_scale(w->A,tf-ti);      
  status=gsl_linalg_exponential_ss(w->A, w->eA, GSL_PREC_DOUBLE);
  assert(status==GSL_SUCCESS);
  for (i=0;i<np;i++){
    b=gsl_matrix_column(Bi,i);
    x=gsl_matrix_column(w->A_B,i);
    status=gsl_linalg_LU_solve(w->LU, w->prm, &b.vector, &x.vector);
    assert(status==GSL_SUCCESS);
  }
  gsl_matrix_memcpy(w->S,Si);
  gsl_matrix_add(w->S,w->A_B);    // S_temp = (Si + A\B)
  gsl_matrix_memcpy(Sf,w->A_B);     // Sf <- A\B
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, w->eA, w->S, -1.0, Sf); // Sf <- exp(At)*(Si + A\B) - A\B
}

/* Peano Baker Series sensitivity approximation */
void PBS
(state_t *initial,
 state_t *final,
 gsl_matrix *PHIf,
 gsl_matrix *PHIb,
 double *work)
{
  assert(initial);
  assert(final);
  int ny=initial->dfdp->size1;
  int np=initial->dfdp->size2;
  gsl_matrix_view temp=gsl_matrix_view_array(work,ny,np);
  double dt=final->t - initial->t;
  transition_matrix_v2
    (initial->dfdy,final->dfdy,
     initial->t,final->t,
     PHIf);
  transition_matrix_v2
    (final->dfdy,initial->dfdy,
     final->t,initial->t,
     PHIb);
  // B = dfdp
  // S = dydp
  // Sf = 1.0*PHI(k+1,k) * (0.5*dt*(1.0*PHI(k,k+1)*B(k+1) + 0.0*B(k+1) + B(k)) + Si ) + 0.0*Sf
  gsl_blas_dgemm
    (CblasNoTrans,CblasNoTrans, 1.0,
     PHIb, final->dfdp,
     0.0,&temp.matrix);
  gsl_matrix_add(&temp.matrix,initial->dfdp);
  gsl_matrix_scale(&temp.matrix,0.5*dt);
  gsl_matrix_add(&temp.matrix,initial->dydp);
  gsl_blas_dgemm
    (CblasNoTrans,CblasNoTrans, 1.0,
     PHIf, &temp.matrix,
     0.0, final->dydp);
}

void update_sensitivity_from_state(gsl_matrix *Sf,state_t *final, int *p, int n){
  int i;
  gsl_vector_view Sf_row, projection_row;
  for (i=0;i<n;i++){
    Sf_row=gsl_matrix_row(Sf,p[i]);
    projection_row=gsl_matrix_row(final->dydp,i);
    gsl_vector_memcpy(&(Sf_row.vector),&(projection_row.vector));
  }
}

/* Calculates the sensitivity matrix using the transition matrix
   PHI. The approximate sensitivity is stored in the `solution`
   structure.*/
void sensitivity_approximation(solution_t *solution)/* solution struct from a numerical solver */{
  assert(solution && solution->Jp && solution->Jy && solution->t);

  int np = solution->Jp->size[0];
  int ny = solution->Jp->size[1];
  int nt = solution->t->size[0];
  int index[3]={0,0,0};
  int j;

  double tf, ti;
  gsl_matrix_view PHIf;
  gsl_matrix_view PHIb;
  gsl_matrix_view Sf,Si,Ji,Jf,Bf,Bi;

  int *p;
  int *n=not_in_steady_state(solution,&p,1e-8,1e-6);
  state_t *initial, *final;
  double *work=malloc(sizeof(double)*ny*np);
  expm_work_t *w=work_mem_alloc(ny,np);
  fprintf(stderr,"[%s] approximating the sensitivity at %i time points\n",__func__,nt);
  for (j=1;j<nt;j++){
    fprintf(stderr,"[%s] %i state variables are not in steady state.\n",__func__,n[j]);
    // create appropriate vector views:
    index[2]=j; // the current time index

    Sf=gsl_matrix_view_array(ndarray_ptr(solution->Sy,index),ny,np);
    Jf=gsl_matrix_view_array(ndarray_ptr(solution->Jy,index),ny,ny);
    Bf=gsl_matrix_view_array(ndarray_ptr(solution->Jp,index),ny,np);
    tf=ndarray_value(solution->t,&(index[2]));
    if (n[j]){
      PHIf=gsl_matrix_view_array(ndarray_ptr(solution->PHIf,index),n[j],n[j]);
      PHIb=gsl_matrix_view_array(ndarray_ptr(solution->PHIb,index),n[j],n[j]);
      
      final=state_projection
	(n[j],&p[ny*j],
	 tf,
	 &Jf.matrix,
	 &Bf.matrix,
	 &Sf.matrix);
    }
    index[2]=j-1; // previous time index
    Si=gsl_matrix_view_array(ndarray_ptr(solution->Sy,index),ny,np);
    Ji=gsl_matrix_view_array(ndarray_ptr(solution->Jy,index),ny,ny);
    Bi=gsl_matrix_view_array(ndarray_ptr(solution->Jp,index),ny,np);
    ti=ndarray_value(solution->t,&(index[2]));
    if (n[j]){
    initial=state_projection
      (n[j],&p[ny*j],
       ti,
       &Ji.matrix,
       &Bi.matrix,
       &Si.matrix);
    }
    /* make a near steady state approximation 
       for the whole state vector */
    near_steady_state_approximation
      (&Sf.matrix,
       &Si.matrix,
       &Ji.matrix,
       &Bi.matrix,
       tf,ti,
       w);
    /* improve on the previous result 
       where steady state condition is not met */
    if(n[j]){
      PBS(initial,final,&PHIf.matrix,&PHIb.matrix,work);
      /* extract the smaller matrix 
         and update the overall result: */
      update_sensitivity_from_state(&Sf.matrix,final,&p[ny*j],n[j]);
      state_projection_free(initial);
      state_projection_free(final);
    }
  }
  free(work);
}
