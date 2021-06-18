#include "solution.h"


solution_t* solution_alloc(int ny, int np, int nt){
  solution_t *solution=malloc(sizeof(solution_t));
  assert(solution);
  int size[3]={np,ny,nt};
  solution->t=ndarray_alloc(1,&size[2]);
  solution->y=ndarray_alloc(2,&size[1]); // list of states
  solution->f=ndarray_alloc(2,&size[1]); // dy/dt

  solution->Sy=ndarray_alloc(3,size);
  solution->Jp=ndarray_alloc(3,size);
  size[0]=ny;
  size[1]=ny;
  size[2]=nt;
  solution->Jy=ndarray_alloc(3,size);
  solution->PHIf=ndarray_alloc(3,size);
  solution->PHIb=ndarray_alloc(3,size);
  solution->status=ndarray_alloc(1,&size[2]);
  return solution;
}

void solution_resize(solution_t *solution, int new_nt){
  assert(solution);
  ndarray_resize(solution->t,new_nt);
  ndarray_resize(solution->y,new_nt);
  ndarray_resize(solution->f,new_nt);
  ndarray_resize(solution->Sy,new_nt);
  ndarray_resize(solution->Jp,new_nt);
  ndarray_resize(solution->Jy,new_nt);
  ndarray_resize(solution->PHIf,new_nt);
  ndarray_resize(solution->PHIb,new_nt);
  ndarray_resize(solution->status,new_nt);
}
