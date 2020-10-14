// The file dbl8_norm_host.h specifies functions to compute
// the 2-norm of a real octo double vector and to normalize a vector.

#ifndef __dbl8_norm_host_h__
#define __dbl8_norm_host_h__

void make_copy
 ( int dim,
   double *orghihihi, double *orglohihi, double *orghilohi, double *orglolohi,
   double *orghihilo, double *orglohilo, double *orghilolo, double *orglololo,
   double *duphihihi, double *duplohihi, double *duphilohi, double *duplolohi,
   double *duphihilo, double *duplohilo, double *duphilolo, double *duplololo);
/*
 * DESCRIPTION :
 *   Makes a copy of the vector in org to the vector dup,
 *   given as eight arrays of the parts of the octo double numbers.
 *
 * ON ENTRY :
 *   dim         dimension of the vectors org and dup;
 *   orghihihi   highest parts of the original vector of dimension dim;
 *   orglohihi   2nd highest parts of the original vector of dimension dim;
 *   orghilohi   3rd highest parts of the original vector of dimension dim;
 *   orglolohi   4th highest parts of the original vector of dimension dim;
 *   orghihilo   4th lowest parts of the original vector of dimension dim;
 *   orglohilo   3rd lowest parts of the original vector of dimension dim;
 *   orghilolo   2nd lowest parts of the original vector of dimension dim;
 *   orglololo   lowest parts of the original vector of dimension dim;
 *   duphihihi   space for the highest parts of a vector of dimension dim;
 *   duplohihi   space for the 2nd highest parts of a vector of dimension dim;
 *   duphilohi   space for the 3rd highest parts of a vector of dimension dim;
 *   duplolohi   space for the 4th highest parts of a vector of dimension dim;
 *   duphihilo   space for the 4th lowest parts of a vector of dimension dim.
 *   duplohilo   space for the 3rd lowest parts of a vector of dimension dim.
 *   duphilolo   space for the 2nd lowest parts of a vector of dimension dim.
 *   duplololo   space for the lowest parts of a vector of dimension dim.
 *
 * ON RETURN :
 *   duphihihi   duplicated vector with same contents as orghihi;
 *   duplohihi   duplicated vector with same contents as orglohi;
 *   duphilohi   duplicated vector with same contents as orghilo;
 *   duplolohi   duplicated vector with same contents as orglolo;
 *   duphihilo   duplicated vector with same contents as orghihi;
 *   duplohilo   duplicated vector with same contents as orglohi;
 *   duphilolo   duplicated vector with same contents as orghilo;
 *   duplololo   duplicated vector with same contents as orglolo. */

void CPU_norm
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the vector in v, given by eight arrays,
 *   of all parts of the octo double numbers.
 *
 * ON ENTRY :
 *   vhihihi     highest parts of a vector of size dim;
 *   vlohihi     second highest parts of a vector of size dim;
 *   vhilohi     third highest parts of a vector of size dim;
 *   vlolohi     fourth highest parts of a vector of size dim;
 *   vhihilo     fourth lowest parts of a vector of size dim;
 *   vlohilo     third lowest parts of a vector of size dim;
 *   vhilolo     second lowest parts of a vector of size dim;
 *   vlololo     lowest parts of a vector of size dim;
 *   dim         dimension of the vector v.
 *
 * ON RETURN :
 *   normhihihi  highest part of the 2-norm of the vector v;
 *   normlohihi  second highest part of the 2-norm of the vector v;
 *   normhilohi  third highest part of the 2-norm of the vector v;
 *   normlolohi  fourth highest part of the 2-norm of the vector v;
 *   normhihilo  fourth lowest part of the 2-norm of the vector v;
 *   normlohilo  third lowest part of the 2-norm of the vector v;
 *   normhilolo  second lowest part of the 2-norm of the vector v;
 *   normlololo  lowest part of the 2-norm of the vector v. */

void CPU_normalize
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, double normhihihi, double normlohihi, double normhilohi,
   double normlolohi, double normhihilo, double normlohilo,
   double normhilolo, double normlololo );
/*
 * DESCRIPTION :
 *   Normalizes the vector in v.
 *
 * ON ENTRY :
 *   vhihihi     highest parts of a real vector, of size dim;
 *   vlohihi     second highest parts of a real vector, of size dim;
 *   vhilohi     third highest parts of a real vector, of size dim;
 *   vlolohi     fourth highest parts of a real vector, of size dim;
 *   vhihilo     fourth lowest parts of a real vector, of size dim;
 *   vlohilo     third lowest parts of a real vector, of size dim;
 *   vhilolo     second lowest parts of a real vector, of size dim;
 *   vlololo     lowest parts of a real vector, of size dim;
 *   dim         dimension of the vector v;
 *   normhihihi  highest part of a nonzero number to divide v with;
 *   normlohihi  second highest  part of a nonzero number to divide v with;
 *   normhilohi  third highest  part of a nonzero number to divide v with;
 *   normlolohi  fourth highest  part of a nonzero number to divide v with;
 *   normhihilo  fourth lowest part of a nonzero number to divide v with.
 *   normlohilo  third lowest part of a nonzero number to divide v with.
 *   normhilolo  second lowest part of a nonzero number to divide v with.
 *   normlololo  lowest part of a nonzero number to divide v with.
 *
 * ON RETURN :
 *   vhihihi     highest parts of the normalized vector;
 *   vlohihi     second highest parts of the normalized vector;
 *   vhilohi     third highest parts of the normalized vector;
 *   vlolohi     fourth highest parts of the normalized vector;
 *   vhihilo     fourth lowest parts of the normalized vector;
 *   vlohilo     third lowest parts of the normalized vector;
 *   vhilolo     second lowest parts of the normalized vector;
 *   vlololo     lowest parts of the normalized vector,
 *               if on entry, norm equals the 2-norm of v, 
 *               then the 2-norm of v on return will be one. */

#endif
