#ifndef SOLUTION_H
#define SOLUTION_H
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "ndarray.h"

typedef struct {
  ndarray *t;
  ndarray *y;
  ndarray *f;
  ndarray *Jy;
  ndarray *Jp; 
  ndarray *Sy;
  ndarray *PHIf;
  ndarray *PHIb;
  ndarray *status;
} solution_t;

solution_t* solution_alloc(int ny, int np, int nt);
void solution_resize(solution_t *solution, int new_nt);
#endif
