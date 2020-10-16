// The file cmplx4_norm_host.h specifies functions to compute
// the 2-norm of a complex vector and to normalize a complex vector,
// in quad double precision.

#ifndef __cmplx4_norm_host_h__
#define __cmplx4_norm_host_h__

void make_copy
 ( int dim,
   double *orgrehihi, double *orgrelohi, double *orgrehilo, double *orgrelolo,
   double *orgimhihi, double *orgimlohi, double *orgimhilo, double *orgimlolo,
   double *duprehihi, double *duprelohi, double *duprehilo, double *duprelolo,
   double *dupimhihi, double *dupimlohi, double *dupimhilo, double *dupimlolo );
/*
 * DESCRIPTION :
 *   Makes a copy of the complex vector with real parts in orgre
 *   (orgrehihi, orgrelohi, orgrehilo, orgrelolo), and imaginary parts 
 *   in orgim (orgimhihi, orgimlohi, orgimhilo, orgimlolo), to the 
 *   complex vector with real parts in dupre (duprehihi, duprelohi, 
 *   duprehilo,duprelolo) and imaginary parts in dupim (dupimhihii,
 *   dupimlohi, dupimhilo, dupimlolo).
 *
 * REQUIRED :
 *   Space has been allocated for all vectors,
 *   to hold at least as many doubles as the value of dim.
 *
 * ON ENTRY :
 *   dim         dimension of all vectors;
 *   orgrehihi   highest real parts of the original vector;
 *   orgrelohi   second highest real parts of the original vector;
 *   orgrehilo   second lowest real parts of the original vector;
 *   orgrelolo   lowest real parts of the original vector;
 *   orgimhihi   highest imaginary parts of the original vector;
 *   orgimlohi   second highest imaginary parts of the original vector;
 *   orgimhilo   second lowest imaginary parts of the original vector;
 *   orgimlolo   lowest imaginary parts of the original vector;
 *   duprehihi   space for highest real parts of a vector;
 *   duprelohi   space for second highest real parts of a vector;
 *   duprehilo   space for second lowest real parts of a vector;
 *   duprelolo   space for lowest real parts of a vector;
 *   dupimhihi   space for highest imaginary parts of a vector;
 *   dupimlohi   space for second highest imaginary parts of a vector;
 *   dupimhilo   space for second lowest imaginary parts of a vector;
 *   dupimlolo   space for lowest imaginary parts of a vector.
 *
 * ON RETURN :
 *   duprehihi   highest real parts of the duplicate vector;
 *   duprelohi   second highest real parts of the duplicate vector;
 *   duprehilo   second lowest real parts of the duplicate vector;
 *   duprelolo   lowest real parts of the duplicate vector;
 *   dupimhihi   highest imaginary parts of the duplicate vector;
 *   dupimlohi   second highest imaginary parts of the duplicate vector;
 *   dupimhilo   second lowest imaginary parts of the duplicate vector;
 *   dupimlolo   lowest imaginary parts of the duplicate vector. */

void CPU_norm
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the complex vector given by its real parts
 *   in (vrehihi, vrelohi, vrehilo, vrelolo) and its imaginary parts
 *   in (vimhhii, vimlohi, vimhilo, vimlolo).
 *
 * REQUIRED : All arrays have size dim.
 *
 * ON ENTRY :
 *   vrehihi   highest real parts of a quad double vector;
 *   vrelohi   second highest real parts of a quad double vector;
 *   vrehilo   second lowest real parts of a quad double vector;
 *   vrelolo   lowest real parts of a quad double vector;
 *   vimhihi   highest imaginary parts of a quad double vector;
 *   vimlohi   second highest imaginary parts of a quad double vector;
 *   vimhilo   second lowest imaginary parts of a quad double vector;
 *   vimlolo   lowest imaginary parts of a quad double vector;
 *   dim       dimension of all arrays.
 *
 * ON RETURN :
 *   normhihi  second highest part of the 2-norm of the complex vector;
 *   normlohi  second highest part of the 2-norm of the complex vector;
 *   normhilo  second lowest part of the 2-norm of the complex vector;
 *   normlolo  lowest part of the 2-norm of the complex vector. */

void CPU_normalize
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim,
   double normhihi, double normlohi, double normhilo, double normlolo );
/*
 * DESCRIPTION :
 *   Normalizes the complex vector given by four arrays for the
 *   real parts (vrehihi, vrelohi, vrehilo, vrelolo) and by four arrays
 *   for the imaginary parts (vimhihi, vimlohi, vimhilo, vimlolo).
 *
 * REQUIRED : All arrays have size dim.
 *
 * ON ENTRY :
 *   vrehihi   highest real parts of a quad double vector;
 *   vrelohi   second highest real parts of a quad double vector;
 *   vrehilo   second lowest real parts of a quad double vector;
 *   vrelolo   lowest real parts of a quad double vector;
 *   vimhihi   highest imaginary parts of a quad double vector;
 *   vimlohi   second highest imaginary parts of a quad double vector;
 *   vimhilo   second lowest imaginary parts of a quad double vector;
 *   vimlo     lowest imaginary parts of a quad double vector;
 *   dim       dimension of the vector v;
 *   normhihi  highest part of the norm to normalize the vector;
 *   normlohi  second highest part of the norm to normalize the vector;
 *   normhilo  second lowest part of the norm to normalize the vector;
 *   normlolo  lowest part of the norm to normalize the vector.
 *
 * ON RETURN :
 *   vrehihi   highest real parts of the normalized vector;
 *   vrelohi   second highest real parts of the normalized vector;
 *   vrehilo   second lowest real parts of the normalized vector;
 *   vrelolo   lowest real parts of the normalized vector;
 *   vimhihi   highest imaginary parts of the normalized vector;
 *   vimlohi   second highest imaginary parts of the normalized vector;
 *   vimhilo   second lowest imaginary parts of the normalized vector;
 *   vimlolo   lowest imaginary parts of the normalized vector;
 *             if on entry, (normhihi, normlohi, normhilo, normlolo) 
 *             equals the 2-norm, then the complex quad double vector
 *             represented by (vrehihi, vrelohi, vrehilo, vrelolo) and
 *             (vimhihi, vimlohi, vimhilo, vimlolo) has 2-norm one. */

#endif
