#include "ndarray.h"

void ndarray_test(){
  int size[3]={4,3,2};
  ndarray *a=ndarray_alloc(3,size);
  int i;
  for (i=0;i<a->prod[2];i++) a->value[i]=i;
  int j[3]={1,1,1};
  double v=ndarray_value(a,j);
  printf("[%s] j=[1,1,1], value=%g (should be %i)\n",__func__,v,17);
  ndarray_print(a,"test array");
  ndarray_resize(a,4);
  ndarray_print(a,"same array resized from 2 to 4 (last index)");
  ndarray_free(a);
  fflush(stdout);
  printf("[%s] ndarray freed.\n",__func__);
}

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


void ndarray_resize(ndarray *a, int new_size){
  assert(new_size>0);
  int r=a->rank;
  a->size[r-1]=new_size;
  if (r>1){  
    a->prod[r-1]=new_size*(a->prod[r-2]);
  } else {
    a->prod[r-1]=new_size;
  }
  a->value=realloc(a->value,sizeof(double)*(a->prod[r-1]));
}

double* ndarray_ptr(ndarray *a, int *index){
  assert(a);
  assert(index);
  int i,j=0,n=1;
  for (i=0;i<a->rank;i++){
    j+=index[i]*n;
    n=a->prod[i];
  }
  assert(j<n);
  return &(a->value[j]);
}

double ndarray_value(ndarray *a, int *index){
  return *(ndarray_ptr(a,index));
}


void ndarray_to_h5(ndarray *a, hid_t loc_id, const char *obj_name){
  herr_t status=0;
  assert(a);
  hsize_t *dims = malloc(sizeof(hsize_t)*a->rank);
  int i;
  int rank=a->rank;
  // hdf5 stores things row-wise, so, it interprets the size differently
  for (i=0;i<rank;i++) dims[i] = (hsize_t) a->size[rank-1-i];
  status=H5LTmake_dataset_double(loc_id, obj_name, (hsize_t) a->rank, dims, a->value);
  assert(status>=0);
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
