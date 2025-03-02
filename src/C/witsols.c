/* The file witsols.c contains the definitions of the functions with
 * prototypes documented in witsols.h. */

#include "witsols.h"
#include <stdio.h>

/* the solvers */

int standard_polysys_solve(int nbtasks, int topdim, int filter, int factor,
                           int verbose) {
  int fail = 0;
  int pars[4];
  double *c;

  pars[0] = nbtasks;
  pars[1] = topdim;
  pars[2] = filter;
  pars[3] = factor;

  fail = _ada_use_c2phc(845, pars, &verbose, c, 0);

  return fail;
}

int standard_laursys_solve(int nbtasks, int topdim, int filter, int factor,
                           int verbose) {
  int fail = 0;
  int pars[4];
  double *c;

  pars[0] = nbtasks;
  pars[1] = topdim;
  pars[2] = filter;
  pars[3] = factor;

  fail = _ada_use_c2phc(846, pars, &verbose, c, 0);

  return fail;
}

int dobldobl_polysys_solve(int nbtasks, int topdim, int filter, int factor,
                           int verbose) {
  int fail = 0;
  int pars[4];
  double *c;

  pars[0] = nbtasks;
  pars[1] = topdim;
  pars[2] = filter;
  pars[3] = factor;

  fail = _ada_use_c2phc(847, pars, &verbose, c, 0);

  return fail;
}

int dobldobl_laursys_solve(int nbtasks, int topdim, int filter, int factor,
                           int verbose) {
  int fail = 0;
  int pars[4];
  double *c;

  pars[0] = nbtasks;
  pars[1] = topdim;
  pars[2] = filter;
  pars[3] = factor;

  fail = _ada_use_c2phc(848, pars, &verbose, c, 0);

  return fail;
}

int quaddobl_polysys_solve(int nbtasks, int topdim, int filter, int factor,
                           int verbose) {
  int fail = 0;
  int pars[4];
  double *c;

  pars[0] = nbtasks;
  pars[1] = topdim;
  pars[2] = filter;
  pars[3] = factor;

  fail = _ada_use_c2phc(849, pars, &verbose, c, 0);

  return fail;
}

int quaddobl_laursys_solve(int nbtasks, int topdim, int filter, int factor,
                           int verbose) {
  int fail = 0;
  int pars[4];
  double *c;

  pars[0] = nbtasks;
  pars[1] = topdim;
  pars[2] = filter;
  pars[3] = factor;

  fail = _ada_use_c2phc(850, pars, &verbose, c, 0);

  return fail;
}

/* extracting solution data */

int copy_standard_polysys_witset(int dim) {
  int *b;
  double *c;

  int fail = _ada_use_c2phc(851, &dim, b, c, 0);

  return fail;
}

int copy_standard_laursys_witset(int dim) {
  int *b;
  double *c;

  int fail = _ada_use_c2phc(852, &dim, b, c, 0);

  return fail;
}

int copy_dobldobl_polysys_witset(int dim) {
  int *b;
  double *c;

  int fail = _ada_use_c2phc(853, &dim, b, c, 0);

  return fail;
}

int copy_dobldobl_laursys_witset(int dim) {
  int *b;
  double *c;

  int fail = _ada_use_c2phc(854, &dim, b, c, 0);

  return fail;
}

int copy_quaddobl_polysys_witset(int dim) {
  int *b;
  double *c;

  int fail = _ada_use_c2phc(855, &dim, b, c, 0);

  return fail;
}

int copy_quaddobl_laursys_witset(int dim) {
  int *b;
  double *c;

  int fail = _ada_use_c2phc(856, &dim, b, c, 0);

  return fail;
}

int clear_standard_witsols(void) {
  int *a, *b;
  double *c;

  int fail = _ada_use_c2phc(857, a, b, c, 0);

  return fail;
}

int clear_dobldobl_witsols(void) {
  int *a, *b;
  double *c;

  int fail = _ada_use_c2phc(858, a, b, c, 0);

  return fail;
}

int clear_quaddobl_witsols(void) {
  int *a, *b;
  double *c;

  int fail = _ada_use_c2phc(859, a, b, c, 0);

  return fail;
}
