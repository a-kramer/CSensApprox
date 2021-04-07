#ifndef NDARRAY_H
#define NDARRAY_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "hdf5.h"
#include "hdf5_hl.h"
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
void ndarray_test();
void ndarray_to_h5(ndarray *a, hid_t loc_id, const char *obj_name);
#endif
