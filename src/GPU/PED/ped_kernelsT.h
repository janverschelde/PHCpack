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
   int ncoefs, char *pos_arr_h_char, char *exp_arr_h_char, complexD<T> *x_h,
   complexD<T> *c_h, complexD<T> *factors_h, complexD<T> *polvalues_h );

#endif
