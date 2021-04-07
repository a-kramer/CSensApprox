#include "solution.h"


solution_t* solution_alloc(int ny, int np, int nt){
  solution_t *solution=malloc(sizeof(solution_t));
  assert(solution);
  solution->t=gsl_vector_alloc(nt);
  solution->y=gsl_matrix_alloc(nt,ny);
  int size[3]={ny,np,nt};
  solution->Sy=ndarray_alloc(3,size);
  size[0]=ny;
  size[1]=ny;
  size[2]=nt;
  solution->Jy=ndarray_alloc(3,size);
  size[0]=np;
  size[1]=ny;
  size[2]=nt;
  solution->Jp=ndarray_alloc(3,size);
  return solution;
}
