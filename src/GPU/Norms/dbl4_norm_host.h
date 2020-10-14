// The file dbl4_norm_host.h specifies functions to compute
// the 2-norm of a real quad double vector and to normalize a vector.

#ifndef __dbl4_norm_host_h__
#define __dbl4_norm_host_h__

void make_copy
 ( int dim,
   double *orghihi, double *orglohi, double *orghilo, double *orglolo,
   double *duphihi, double *duplohi, double *duphilo, double *duplolo);
/*
 * DESCRIPTION :
 *   Makes a copy of the vector in org to the vector dup,
 *   given as vectors of highest, second highest, second lowest,
 *   and lowest parts.
 *
 * ON ENTRY :
 *   dim       dimension of the vectors org and dup;
 *   orghihi   highest parts of the original vector of dimension dim;
 *   orglohi   2nd highest parts of the original vector of dimension dim;
 *   orghilo   2nd lowest parts of the original vector of dimension dim;
 *   orglolo   lowest parts of the original vector of dimension dim;
 *   duphihi   space for the highest parts of a vector of dimension dim;
 *   duplohi   space for the 2nd highest parts of a vector of dimension dim;
 *   duphilo   space for the 2nd lowest parts of a vector of dimension dim.
 *   duplolo   space for the lowest parts of a vector of dimension dim.
 *
 * ON RETURN :
 *   duphihi   duplicated vector with same contents as orghihi;
 *   duplohi   duplicated vector with same contents as orglohi;
 *   duphilo   duplicated vector with same contents as orghilo.
 *   duplolo   duplicated vector with same contents as orglolo. */

void CPU_norm
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo, int dim,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the vector in v, given by four arrays,
 *   of highest, second highest, second lowest and lowest parts
 *   of quad double numbers.
 *
 * ON ENTRY :
 *   vhihi     highest parts of a vector of size dim;
 *   vlohi     second highest parts of a vector of size dim;
 *   vhilo     second lowest parts of a vector of size dim;
 *   vlolo     lowest parts of a vector of size dim;
 *   dim       dimension of the vector v.
 *
 * ON RETURN :
 *   normhihi  highest part of the 2-norm of the vector v;
 *   normlohi  second highest part of the 2-norm of the vector v;
 *   normhilo  second lowest part of the 2-norm of the vector v;
 *   normlolo  lowest part of the 2-norm of the vector v. */

void CPU_normalize
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo, int dim,
   double normhihi, double normlohi, double normhilo, double normlolo );
/*
 * DESCRIPTION :
 *   Normalizes the vector in v.
 *
 * ON ENTRY :
 *   vhihi     highest parts of a real vector, of size dim;
 *   vlohi     second highest parts of a real vector, of size dim;
 *   vhilo     second lowest parts of a real vector, of size dim;
 *   vlolo     lowest parts of a real vector, of size dim;
 *   dim       dimension of the vector v;
 *   normhihi  highest part of a nonzero number to divide v with;
 *   normlohi  second highest  part of a nonzero number to divide v with;
 *   normhilo  second lowest part of a nonzero number to divide v with.
 *   normlolo  lowest part of a nonzero number to divide v with.
 *
 * ON RETURN :
 *   vhihi     highest parts of the normalized vector;
 *   vlohi     second highest parts of the normalized vector;
 *   vhilo     second lowest parts of the normalized vector;
 *   vlolo     lowest parts of the normalized vector,
 *             if on entry, norm equals the 2-norm of v, 
 *             then the 2-norm of v on return will be one. */

#endif
