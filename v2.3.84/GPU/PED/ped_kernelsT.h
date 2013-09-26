// Kernels to evaluate and differentiate polynomials in several variables.
// Originally written by Genady Yoffe with modifications by Jan Verschelde.

#ifndef __KERNELS__
#define __KERNELS__

#include <iostream>
#include <gqd_type.h>
#include "DefineType.h"


template<class T>
void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int mode, int m,
   int ncoefs, char *pos_arr_h_char, char *exp_arr_h_char, complex<T> *x_h,
   complex<T> *c_h, complex<T> *factors_h, complex<T> *polvalues_h );

#endif
