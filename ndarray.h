#ifndef NDARRAY_H
#define NDARRAY_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* these words are the same as in the GNU scientific library */
typedef struct {
  int rank;
  int *size;
  int *prod;
  double *value;
} ndarray;

ndarray* ndarray_alloc(int rank, int *size);
double ndarray_value(ndarray *a, int *index);
double* ndarray_ptr(ndarray *a, int *index);
void ndarray_print(ndarray *a, char *label);
void ndarray_free(ndarray *a);

#endif
