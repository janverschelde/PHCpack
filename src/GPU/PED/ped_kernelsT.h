// Kernels to evaluate and differentiate polynomials in several variables.
// Originally written by Genady Yoffe with modifications by Jan Verschelde.

#ifndef __KERNELS__
#define __KERNELS__

#include <iostream>
#include <gqd_type.h>
#include "DefineType.h"

template<class realD>
void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m,
   int ncoefs, char *pos_arr_h_char, char *exp_arr_h_char,
   complexD<realD> *x_h, complexD<realD> *c_h,
   complexD<realD> *factors_h, complexD<realD> *polvalues_h );

#endif
