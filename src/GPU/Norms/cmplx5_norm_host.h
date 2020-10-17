// The file cmplx5_norm_host.h specifies functions to compute
// the 2-norm of a complex vector and to normalize a complex vector,
// in penta double precision.

#ifndef __cmplx5_norm_host_h__
#define __cmplx5_norm_host_h__

void make_copy
 ( int dim,
   double *orgretb, double *orgreix, double *orgremi, double *orgrerg,
   double *orgrepk,
   double *orgimtb, double *orgimix, double *orgimmi, double *orgimrg,
   double *orgimpk,
   double *dupretb, double *dupreix, double *dupremi, double *duprerg,
   double *duprepk,
   double *dupimtb, double *dupimix, double *dupimmi, double *dupimrg,
   double *dupimpk );
/*
 * DESCRIPTION :
 *   Makes a copy of the complex vector with real parts in orgre
 *   (orgretb, orgreix, orgremi, orgrerg, orgrepk), and imaginary parts 
 *   in orgim (orgimtb, orgimix, orgimmi, orgimrg, orgimpk), to the 
 *   complex vector with real parts in dupre (dupretb, dupreix, 
 *   dupremi, duprerg, duprepk) and imaginary parts in dupim (dupimtbi,
 *   dupimix, dupimmi, dupimrg, dupimpk).
 *
 * REQUIRED :
 *   Space has been allocated for all vectors,
 *   to hold at least as many doubles as the value of dim.
 *
 * ON ENTRY :
 *   dim       dimension of all vectors;
 *   orgretb   highest real parts of the original vector;
 *   orgreix   second highest real parts of the original vector;
 *   orgremi   middle real parts of the original vector;
 *   orgrerg   second lowest real parts of the original vector;
 *   orgrepk   lowest real parts of the original vector;
 *   orgimtb   highest imaginary parts of the original vector;
 *   orgimix   second highest imaginary parts of the original vector;
 *   orgimmi   middle imaginary parts of the original vector;
 *   orgimrg   second lowest imaginary parts of the original vector;
 *   orgimpk   lowest imaginary parts of the original vector;
 *   dupretb   space for highest real parts of a vector;
 *   dupreix   space for second highest real parts of a vector;
 *   dupremi   space for middle real parts of a vector;
 *   duprerg   space for second lowest real parts of a vector;
 *   duprepk   space for lowest real parts of a vector;
 *   dupimtb   space for highest imaginary parts of a vector;
 *   dupimix   space for second highest imaginary parts of a vector;
 *   dupimmi   space for middle imaginary parts of a vector;
 *   dupimrg   space for second lowest imaginary parts of a vector;
 *   dupimpk   space for lowest imaginary parts of a vector.
 *
 * ON RETURN :
 *   dupretb   highest real parts of the duplicate vector;
 *   dupreix   second highest real parts of the duplicate vector;
 *   dupremi   middle real parts of the duplicate vector;
 *   duprerg   second lowest real parts of the duplicate vector;
 *   duprepk   lowest real parts of the duplicate vector;
 *   dupimtb   highest imaginary parts of the duplicate vector;
 *   dupimix   second highest imaginary parts of the duplicate vector;
 *   dupimmi   middle imaginary parts of the duplicate vector;
 *   dupimrg   second lowest imaginary parts of the duplicate vector;
 *   dupimpk   lowest imaginary parts of the duplicate vector. */

void CPU_norm
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the complex vector given by its real parts
 *   in (vretb, vreix, vremi, vrerg, vrepk) and its imaginary parts
 *   in (vimhhii, vimix, vimmi, vimrg, vimpk).
 *
 * REQUIRED : All arrays have size dim.
 *
 * ON ENTRY :
 *   vretb    highest real parts of a penta double vector;
 *   vreix    second highest real parts of a penta double vector;
 *   vremi    middle real parts of a penta double vector;
 *   vrerg    second lowest real parts of a penta double vector;
 *   vrepk    lowest real parts of a penta double vector;
 *   vimtb    highest imaginary parts of a penta double vector;
 *   vimix    second highest imaginary parts of a penta double vector;
 *   vimmi    middle imaginary parts of a penta double vector;
 *   vimrg    second lowest imaginary parts of a penta double vector;
 *   vimpk    lowest imaginary parts of a penta double vector;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   normtb   highest part of the 2-norm of the complex vector;
 *   normix   second highest part of the 2-norm of the complex vector;
 *   normmi   middle part of the 2-norm of the complex vector;
 *   normrg   second lowest part of the 2-norm of the complex vector;
 *   normpk   lowest part of the 2-norm of the complex vector. */

void CPU_normalize
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim,
   double normtb, double normix, double normmi, double normrg, double normpk );
/*
 * DESCRIPTION :
 *   Normalizes the complex vector given by five arrays for the
 *   real parts (vretb, vreix, vremi, vrerg, vrepk) and by five arrays
 *   for the imaginary parts (vimtb, vimix, vimmi, vimrg, vimpk).
 *
 * REQUIRED : All arrays have size dim.
 *
 * ON ENTRY :
 *   vretb    highest real parts of a penta double vector;
 *   vreix    second highest real parts of a penta double vector;
 *   vremi    middle real parts of a penta double vector;
 *   vrerg    second lowest real parts of a penta double vector;
 *   vrepk    lowest real parts of a penta double vector;
 *   vimtb    highest imaginary parts of a penta double vector;
 *   vimix    second highest imaginary parts of a penta double vector;
 *   vimmi    middle imaginary parts of a penta double vector;
 *   vimrg    second lowest imaginary parts of a penta double vector;
 *   vimpk    lowest imaginary parts of a penta double vector;
 *   dim      dimension of the vector v;
 *   normtb   highest part of the norm to normalize the vector;
 *   normix   second highest part of the norm to normalize the vector;
 *   normmi   middle part of the norm to normalize the vector;
 *   normrg   second lowest part of the norm to normalize the vector.
 *   normpk   lowest part of the norm to normalize the vector.
 *
 * ON RETURN :
 *   vretb    highest real parts of the normalized vector;
 *   vreix    second highest real parts of the normalized vector;
 *   vremi    middle real parts of the normalized vector;
 *   vrerg    second lowest real parts of the normalized vector;
 *   vrepk    lowest real parts of the normalized vector;
 *   vimtb    highest imaginary parts of the normalized vector;
 *   vimix    second highest imaginary parts of the normalized vector;
 *   vimmi    middle imaginary parts of the normalized vector;
 *   vimrg    second lowest imaginary parts of the normalized vector;
 *   vimpk    lowest imaginary parts of the normalized vector;
 *            if on entry, (normtb, normix, normmi, normrg, normpk) 
 *            equals the 2-norm, then the complex penta double vector
 *            represented by (vretb, vreix, vremi, vrerg, vrepk) and
 *            (vimtb, vimix, vimmi, vimrg, vimpk) has 2-norm one. */

#endif
