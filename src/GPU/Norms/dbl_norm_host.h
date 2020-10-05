// The file dbl_norm_host.h specifies functions to compute
// the 2-norm of a real double vector and to normalize a vector.

#ifndef __dbl_norm_host_h__
#define __dbl_norm_host_h__

void make_copy ( int dim, double* org, double* dup );
/*
 * DESCRIPTION :
 *   Makes a copy of the vector in org to the vector dup.
 *
 * ON ENTRY :
 *   dim       dimension of the vectors org and dup;
 *   org       original vector of dimension dim;
 *   dup       allocated space for a vector of dimension dim.
 *
 * ON RETURN :
 *   dup       duplicated vector with same contents as org. */

void CPU_norm ( double* v, int dim, double* twonorm );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the vector in v.
 *
 * ON ENTRY :
 *   v         real vector of size dim;
 *   dim       dimension of the vector v.
 *
 * ON RETURN :
 *   twonorm   the 2-norm of the vector v. */

void CPU_normalize ( double* v, int dim, double norm );
/*
 * DESCRIPTION :
 *   Normalizes the vector in v.
 *
 * ON ENTRY :
 *   v         vector of doubles, of size dim;
 *   dim       dimension of the vector v;
 *   norm      nonzero number to divide every element in v with.
 *
 * ON RETURN :
 *   v         vector with every element divided by norm,
 *             if on entry, norm equals the 2-norm of v, 
 *             then the 2-norm of v on return will be one. */

#endif
