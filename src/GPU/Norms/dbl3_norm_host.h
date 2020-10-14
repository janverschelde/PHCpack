// The file dbl3_norm_host.h specifies functions to compute
// the 2-norm of a real triple double vector and to normalize a vector.

#ifndef __dbl3_norm_host_h__
#define __dbl3_norm_host_h__

void make_copy
 ( int dim, double *orghi, double *orgmi, double *orglo,
            double *duphi, double *dupmi, double *duplo );
/*
 * DESCRIPTION :
 *   Makes a copy of the vector in org to the vector dup,
 *   given as vectors of high, middle, and low parts.
 *
 * ON ENTRY :
 *   dim       dimension of the vectors org and dup;
 *   orghi     high parts of the original vector of dimension dim;
 *   orgmi     middle parts of the original vector of dimension dim;
 *   orglo     low parts of the original vector of dimension dim;
 *   duphi     space for the high parts of a vector of dimension dim;
 *   dupmi     space for the middle parts of a vector of dimension dim;
 *   duplo     space for the low parts of a vector of dimension dim.
 *
 * ON RETURN :
 *   duphi     duplicated vector with same contents as orghi;
 *   dupmi     duplicated vector with same contents as orgmi;
 *   duplo     duplicated vector with same contents as orglo. */

void CPU_norm
 ( double *vhi, double *vmi, double *vlo, int dim,
   double *normhi, double *normmi, double *normlo );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the vector in v, given by three vectors,
 *   of high, middle and low parts of triple doubles.
 *
 * ON ENTRY :
 *   vhi       high parts of a vector of size dim;
 *   vmi       middle parts of a vector of size dim;
 *   vlo       low parts of a vector of size dim;
 *   dim       dimension of the vector v.
 *
 * ON RETURN :
 *   normhi    high part of the 2-norm of the vector v;
 *   normmi    middle part of the 2-norm of the vector v;
 *   normlo    low part of the 2-norm of the vector v. */

void CPU_normalize
 ( double *vhi, double *vmi, double *vlo, int dim,
   double normhi, double normmi, double normlo );
/*
 * DESCRIPTION :
 *   Normalizes the vector in v.
 *
 * ON ENTRY :
 *   vhi       high parts of a real vector, of size dim;
 *   vmi       middle parts of a real vector, of size dim;
 *   vlo       low parts of a real vector, of size dim;
 *   dim       dimension of the vector v;
 *   normhi    high part of a nonzero number to divide v with;
 *   normmi    middle part of a nonzero number to divide v with;
 *   normlo    low part of a nonzero number to divide v with.
 *
 * ON RETURN :
 *   vhi       high parts of the normalized vector;
 *   vmi       middle parts of the normalized vector;
 *   vlo       low parts of the normalized vector,
 *             if on entry, norm equals the 2-norm of v, 
 *             then the 2-norm of v on return will be one. */

#endif
