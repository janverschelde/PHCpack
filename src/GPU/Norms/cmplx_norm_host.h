// The file cmplx_norm_host.h specifies functions to compute
// the 2-norm of a complex vector and to normalize a complex vector.

#ifndef __cmplx_norm_host_h__
#define __cmplx_norm_host_h__

void make_copy
 ( int dim, double* orgre, double* orgim, double* dupre, double* dupim );
/*
 * DESCRIPTION :
 *   Makes a copy of the complex vector with real parts in orgre
 *   and imaginary parts in orgim to the complex vector with real
 *   parts in dupre and imaginary parts in dupim.
 *
 * REQUIRED :
 *   Space has been allocated for all vectors,
 *   to hold at least as many doubles as the value of dim.
 *
 * ON ENTRY :
 *   dim       dimension of all vectors;
 *   orgre     real parts of the original vector of dimension dim;
 *   orgim     imaginary parts of the original vector of dimension dim;
 *   dupre     space for real parts of a vector of dimension dim;
 *   dupim     space for imaginary parts of a vector of dimension dim.
 *
 * ON RETURN :
 *   dupre     real parts of the vector with same contents as org;
 *   dupim     imaginary parts of the vector with same contents as org. */

void CPU_norm ( double* vre, double* vim, int dim, double* twonorm );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the complex vector given as a vector
 *   of real parts vre and a vector of imaginary parts in vim.
 *
 * ON ENTRY :
 *   vre       real parts of a vector of complex numbers of size dim;
 *   vim       imaginary parts of a vector of complex numbers of size dim;
 *   dim       dimension of the vector v.
 *
 * ON RETURN :
 *   twonorm   the 2-norm of the complex vector given by vre and vim. */

void CPU_normalize ( double* vre, double* vim, int dim, double norm );
/*
 * DESCRIPTION :
 *   Normalizes the complex vector given as a vector of real parts vre
 *   and a vector of imaginary parts vim.
 *
 * ON ENTRY :
 *   vre       real parts of a vector of doubles, of size dim;
 *   vim       imaginary parts of a vector of doubles, of size dim;
 *   dim       dimension of the vector v;
 *   norm      nonzero number to divide every element in v with.
 *
 * ON RETURN :
 *   vre       real parts of the vector with every element divided by norm;
 *   vim       imaginary parts of the normalized vector,
 *             if on entry, norm equals the 2-norm of the vector,
 *             then the 2-norm of v on return will be one. */

#endif
