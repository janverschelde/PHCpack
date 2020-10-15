// The file cmplx3_norm_host.h specifies functions to compute
// the 2-norm of a complex vector and to normalize a complex vector,
// in triple double precision.

#ifndef __cmplx3_norm_host_h__
#define __cmplx3_norm_host_h__

void make_copy
 ( int dim,
   double *orgrehi, double *orgremi, double *orgrelo,
   double *orgimhi, double *orgimmi, double *orgimlo,
   double *duprehi, double *dupremi, double *duprelo,
   double *dupimhi, double *dupimmi, double *dupimlo );
/*
 * DESCRIPTION :
 *   Makes a copy of the complex vector with real parts
 *   in orgre (orgrehi, orgremi, orgrelo), and imaginary parts
 *   in orgim (orgimhi, orgimmi, orgimlo), to the complex vector
 *   with real parts in dupre (duprehi, dupremi, duprelo), and
 *   imaginary parts in dupim (dupimhi, dupimmi, dupimlo).
 *
 * REQUIRED :
 *   Space has been allocated for all vectors,
 *   to hold at least as many doubles as the value of dim.
 *
 * ON ENTRY :
 *   dim       dimension of all vectors;
 *   orgrehi   high real parts of the original vector;
 *   orgremi   middle real parts of the original vector;
 *   orgrelo   low real parts of the original vector;
 *   orgimhi   high imaginary parts of the original vector;
 *   orgimmi   middle imaginary parts of the original vector;
 *   orgimlo   low imaginary parts of the original vector;
 *   duprehi   space for high real parts of a vector;
 *   dupremi   space for middle real parts of a vector;
 *   duprelo   space for low real parts of a vector;
 *   dupimhi   space for high imaginary parts of a vector;
 *   dupimmi   space for middle imaginary parts of a vector;
 *   dupimlo   space for low imaginary parts of a vector.
 *
 * ON RETURN :
 *   duprehi   high real parts of the duplicate vector to org;
 *   dupremi   middle real parts of the duplicate vector to org;
 *   duprelo   low real parts of the duplicate vector to org;
 *   dupimhi   high imaginary parts of the duplicate vector to org;
 *   dupimmi   middle imaginary parts of the duplicate vector to org;
 *   dupimlo   low imaginary parts of the duplicate vector to org. */

void CPU_norm
 ( double *vrehi, double *vremi, double *vrelo,
   double *vimhi, double *vimmi, double *vimlo, int dim,
   double *normhi, double *normmi, double *normlo );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the complex vector given by its real parts in
 *   (vrehi, vremi, vrelo) and its imaginary parts in (vimhi, vimmi, vimlo).
 *
 * ON ENTRY :
 *   vrehi     high real parts of a triple double vector, of size dim;
 *   vremi     middle real parts of a triple double vector, of size dim;
 *   vrelo     low real parts of a triple double vector, of size dim;
 *   vremi     middle real parts of a triple double vector, of size dim;
 *   vimhi     high imaginary parts of a triple double vector, of size dim;
 *   vimlo     low imaginary parts of a triple double vector, of size dim;
 *   dim       dimension of all arrays.
 *
 * ON RETURN :
 *   normhi    high part of the 2-norm of the complex vector;
 *   normmi    middle part of the 2-norm of the complex vector;
 *   normlo    low part of the 2-norm of the complex vector. */

void CPU_normalize
 ( double *vrehi, double *vremi, double *vrelo,
   double *vimhi, double *vimmi, double *vimlo, int dim,
   double normhi, double normmi, double normlo );
/*
 * DESCRIPTION :
 *   Normalizes the complex vector given by three arrays of 
 *   real parts vre (vrehi, vremi, vrelo) and three arrays of 
 *   imaginary parts vim (vimhi, vimmi, vimlo).
 *
 * ON ENTRY :
 *   vrehi     high real parts of a triple double vector, of size dim;
 *   vremi     middle real parts of a triple double vector, of size dim;
 *   vrelo     low real parts of a triple double vector, of size dim;
 *   vimhi     high imaginary parts of a triple double vector, of size dim;
 *   vimmi     middle imaginary parts of a triple double vector, of size dim;
 *   vimlo     low imaginary parts of a triple double vector, of size dim;
 *   dim       dimension of the vector v;
 *   normhi    high part of the norm to normalize the vector;
 *   normmi    high part of the norm to normalize the vector;
 *   normhi    low part of the norm to normalize the vector.
 *
 * ON RETURN :
 *   vrehi     high real parts of the normalized vector;
 *   vremi     middle real parts of the normalized vector;
 *   vrelo     low real parts of the normalized vector;
 *   vimhi     high imaginary parts of the normalized vector;
 *   vimmi     middle imaginary parts of the normalized vector;
 *   vimlo     low imaginary parts of the normalized vector;
 *             if on entry, (normhi, normmi, normlo) equals the 2-norm,
 *             then the complex triple double vector represented by
 *             (vrehi, vremi, vrelo) and (vimhi, vimimi, vimlo)
 *             has 2-norm one. */

#endif
