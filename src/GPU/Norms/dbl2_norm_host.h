// The file dbl2_norm_host.h specifies functions to compute
// the 2-norm of a real double double vector and to normalize a vector.

#ifndef __dbl2_norm_host_h__
#define __dbl2_norm_host_h__

void make_copy
 ( int dim, double* orghi, double* orglo, double* duphi, double* duplo );
/*
 * DESCRIPTION :
 *   Makes a copy of the vector in org to the vector dup,
 *   given as vectors of high and low parts.
 *
 * ON ENTRY :
 *   dim       dimension of the vectors org and dup;
 *   orghi     high parts of the original vector of dimension dim;
 *   orglo     low parts of the original vector of dimension dim;
 *   duphi     space for the high parts of a vector of dimension dim;
 *   duplo     space for the low parts of a vector of dimension dim.
 *
 * ON RETURN :
 *   duphi     duplicated vector with same contents as orghi;
 *   duplo     duplicated vector with same contents as orglo. */

void CPU_norm
 ( double* vhi, double* vlo, int dim, double* normhi, double* normlo );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the vector in v,
 *   given as a vector of high parts and a vector of low parts.
 *
 * ON ENTRY :
 *   vhi       high parts of a vector of size dim;
 *   vlo       low parts of a vector of size dim;
 *   dim       dimension of the vector v.
 *
 * ON RETURN :
 *   normhi    high part of the 2-norm of the vector v;
 *   normlo    high part of the 2-norm of the vector v. */

void CPU_normalize
 ( double* vhi, double* vlo, int dim, double normhi, double normlo );
/*
 * DESCRIPTION :
 *   Normalizes the vector in v.
 *
 * ON ENTRY :
 *   vhi       high parts of a real vector, of size dim;
 *   vlo       high parts of a real vector, of size dim;
 *   dim       dimension of the vector v;
 *   normhi    high part of a nonzero number to divide v with;
 *   normlo    low part of a nonzero number to divide v with.
 *
 * ON RETURN :
 *   vhi       high parts of the normalized vector;
 *   vlo       low parts of the normalized vector,
 *             if on entry, norm equals the 2-norm of v, 
 *             then the 2-norm of v on return will be one. */

#endif
