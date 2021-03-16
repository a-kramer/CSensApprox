#include "ndarray.h"


ndarray* ndarray_alloc(int rank, int *size){
  assert(rank>0);
  assert(size);
  ndarray *A=malloc(sizeof(ndarray));
  A->rank=rank;
  A->size=malloc(sizeof(int)*rank);
  A->prod=malloc(sizeof(int)*rank);
  
  memcpy(A->size,size,sizeof(int)*rank);
  double p=1;
  int i;
  for (i=0;i<rank;i++) {
    p*=size[i];
    A->prod[i]=p;
  }
  A->value=calloc(p,sizeof(double));
  return A;
}

double ndarray_value(ndarray *a, int *index){
  assert(a);
  assert(index);
  int rank=a->rank;
  int *size=a->size;
  int i,j=0,n=1;
  for (i=0;i<rank;i++){
    j+=size[i]*n;
    n=a->prod[i];
  }
  return a->value[j];
}

double* ndarray_ptr(ndarray *a, int *index){
  assert(a);
  assert(index);
  int rank=a->rank;
  int *size=a->size;
  int i,j=0,n=1;
  for (i=0;i<rank;i++){
    j+=size[i]*n;
    n=a->prod[i];
  }
  return &(a->value[j]);
}


void ndarray_print(ndarray *a, char *label){
  assert(a);
  if (label) printf("#[%s] %s\n",__func__,label);
  else label="value";
  int rank=a->rank;
  int *size=a->size;
  int i,j,k,I,J,K;
  I=size[0];
  J=rank>1?size[1]:1;
  K=rank>2?a->prod[rank-1]/(I*J):1;
  if (rank>2) printf("# all indices beyond i and j will be summarized into one (k)\n");
  for (k=0;k<K;k++){
    printf("%s(:,:,%i) = \n",label,k);
    for (i=0;i<I;i++){
      for (j=0;j<J;j++){
	printf("%g ",a->value[i+j*I+k*I*J]);
      }
      printf("\n");
    }
    printf("\n");
  }
}

void ndarray_free(ndarray *a){
  if (a){
    free(a->size);
    free(a->prod);
    free(a->value);
    free(a);
  }
}
