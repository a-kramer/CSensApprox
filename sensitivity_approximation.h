#ifndef SENS_APPROX_H
#define SENS_APPROX_H
#include <stdlib.h>
#include <stdio.h>
#include "solution.h"

typedef struct {
  double t;
  gsl_matrix *dfdy; // A in the paper
  gsl_matrix *dfdp; // B in the paper
  gsl_matrix *dydp; // sensitivity
  gsl_matrix *PHI_fwd; // transition matrix forward
  gsl_matrix *PHI_bwd; // transition matrix backward
} state_t;

typedef struct {
  gsl_matrix *A;
  gsl_matrix *eA;
  gsl_matrix *LU;
  gsl_permutation *prm;
  int sign;
  gsl_matrix *A_B;
  gsl_matrix *S;
} expm_work_t;

int* is_in_steady_state(solution_t *solution,double abs_tol, double rel_tol);
void sensitivity_approximation(solution_t *solution);
#endif