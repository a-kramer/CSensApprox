#ifndef SOLUTION_H
#define SOLUTION_H
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "ndarray.h"

typedef struct {
  gsl_vector *t;
  gsl_matrix *y;
  ndarray *Jy;
  ndarray *Jp; 
  ndarray *Sy;
} solution_t;

solution_t* solution_alloc(int ny, int np, int nt);
#endif
