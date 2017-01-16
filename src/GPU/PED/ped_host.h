// prototype for the polynomial evaluation and differention on the host

#ifndef _PED_HOST_
#define _PED_HOST_

#include <qd/qd_real.h>
#include "complexD.h"
#include "complexH.h"

void CPU_evaldiff
 ( int dim, int NM, int NV, int deg, int r, int m, int *pos, int *exp,
   complexH<double> *c, complexH<double> *x,
   complexH<double> *factors_h, complexH<double> **polvalues_h );
void CPU_evaldiff
 ( int dim, int NM, int NV, int deg, int r, int m, int *pos, int *exp,
   complexH<dd_real> *c, complexH<dd_real> *x,
   complexH<dd_real> *factors_h, complexH<dd_real> **polvalues_h );
void CPU_evaldiff
 ( int dim, int NM, int NV, int deg, int r, int m, int *pos, int *exp,
   complexH<qd_real> *c, complexH<qd_real> *x,
   complexH<qd_real> *factors_h, complexH<qd_real> **polvalues_h );
/*
 * The CPU is used to evaluate a system and its Jacobian matrix,
 * in double, double double, and quad double precision.
 *
 * ON ENTRY :
 *   dim          dimension of the problem;
 *   NM           number of monomials;
 *   NV           number of variables;
 *   deg          highest degree of the variables;
 *   r            frequency of the runs,
 *   m            NM divided by dim
 *   pos          indicate the participating variables in the monomials;
 *   exp          the exponents of the participating variables; 
 *   c            coefficients of the monomials;
 *   x            point where to evaluate.
 *
 * ON RETURN :
 *   factors_h    the factors common to evaluation and differentiation;
 *   polvalues_h  are the values of polynomials and their derivatives. */

#endif
