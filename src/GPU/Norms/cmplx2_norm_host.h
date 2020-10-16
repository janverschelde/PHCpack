// The file cmplx2_norm_host.h specifies functions to compute
// the 2-norm of a complex vector and to normalize a complex vector,
// in double double precision.

#ifndef __cmplx2_norm_host_h__
#define __cmplx2_norm_host_h__

void make_copy
 ( int dim,
   double *orgrehi, double *orgrelo, double *orgimhi, double *orgimlo,
   double *duprehi, double *duprelo, double *dupimhi, double *dupimlo );
/*
 * DESCRIPTION :
 *   Makes a copy of the complex vector with real parts in orgre
 *   (orgrehi, orgrelo), and imaginary parts in orgim (orgimhi, orgimlo),
 *   to the complex vector with real parts in dupre (duprehi, duprelo)
 *   and imaginary parts in dupim (dupimhi, dupimlo).
 *
 * REQUIRED :
 *   Space has been allocated for all vectors,
 *   to hold at least as many doubles as the value of dim.
 *
 * ON ENTRY :
 *   dim       dimension of all vectors;
 *   orgrehi   high real parts of the original vector;
 *   orgrelo   low real parts of the original vector;
 *   orgimhi   high imaginary parts of the original vector;
 *   orgimlo   low imaginary parts of the original vector;
 *   duprehi   space for high real parts of a vector;
 *   duprelo   space for low real parts of a vector;
 *   dupimhi   space for high imaginary parts of a vector.
 *   dupimlo   space for low imaginary parts of a vector.
 *
 * ON RETURN :
 *   duprehi   high real parts of the duplicate vector to org;
 *   duprelo   low real parts of the duplicate vector to org;
 *   dupimhi   high imaginary parts of the duplicate vector to org;
 *   dupimlo   low imaginary parts of the duplicate vector to org. */

void CPU_norm
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   int dim, double *normhi, double *normlo );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the complex vector given by its real parts
 *   in (vrehi, vrelo) and its imaginary parts in (vimhi, vimlo).
 *
 * ON ENTRY :
 *   vrehi     high real parts of a double double vector, of size dim;
 *   vrelo     low real parts of a double double vector, of size dim;
 *   vimhi     high imaginary parts of a double double vector, of size dim;
 *   vimlo     low imaginary parts of a double double vector, of size dim;
 *   dim       dimension of all arrays.
 *
 * ON RETURN :
 *   normhi    high part of the 2-norm of the complex vector;
 *   normlo    low part of the 2-norm of the complex vector. */

void CPU_normalize
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   int dim, double normhi, double normlo );
/*
 * DESCRIPTION :
 *   Normalizes the complex vector given as a vector of real parts 
 *   (vrehi, vrelo) and a vector of imaginary parts (vimhi, vimlo).
 *
 * ON ENTRY :
 *   vrehi     high real parts of a double double vector, of size dim;
 *   vrelo     low real parts of a double double vector, of size dim;
 *   vimhi     high imaginary parts of a double double vector, of size dim;
 *   vimlo     low imaginary parts of a double double vector, of size dim;
 *   dim       dimension of the vector v;
 *   normhi    high part of the norm to normalize the vector;
 *   normlo    low part of the norm to normalize the vector.
 *
 * ON RETURN :
 *   vrehi     high real parts of the normalized vector;
 *   vrelo     low real parts of the normalized vector;
 *   vimhi     high imaginary parts of the normalized vector;
 *   vimlo     low imaginary parts of the normalized vector;
 *             if on entry, (normhi, normlo) equals the 2-norm,
 *             then the complex double double vector represented by
 *             (vrehi, vrelo) and (vimhi, vimlo) has 2-norm one. */

#endif
