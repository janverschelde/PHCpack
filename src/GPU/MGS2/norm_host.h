// defines code for execution on the host

#ifndef _norm_h
#define _norm_h

#include "DefineType.h"
#include "complex.h"
#include "complexH.h"

void make_copy ( int dim, complex<T>* org, complex<T>* dup );
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

void make_copyH ( int dim, complexH<T1>* org, complexH<T1>* dup );
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

void CPU_norm ( complexH<T1>* v, int dim, T1* twonorm );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the complex vector in v.
 *
 * ON ENTRY :
 *   v         vector of complex numbers of size dim;
 *   dim       dimension of the vector v.
 *
 * ON RETURN :
 *   twonorm   the 2-norm of the vector v. */

void CPU_normalize ( complexH<T1>* v, int dim, T1 norm );
/*
 * DESCRIPTION :
 *   Normalizes the complex vector in v.
 *
 * ON ENTRY :
 *   v         vector of complex numbers of size dim;
 *   dim       dimension of the vector v;
 *   norm      nonzero number to divide every element in v with.
 *
 * ON RETURN :
 *   v         vector with every element divided by norm,
 *             if on entry, norm equals the 2-norm of v, 
 *             then the 2-norm of v on return will be one. */

#endif
