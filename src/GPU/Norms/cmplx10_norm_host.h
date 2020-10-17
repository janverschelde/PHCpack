// The file cmplx10_norm_host.h specifies functions to compute
// the 2-norm of a complex vector and to normalize a complex vector,
// in deca double precision.

#ifndef __cmplx10_norm_host_h__
#define __cmplx10_norm_host_h__

void make_copy
 ( int dim,
   double *orgrertb, double *orgrerix, double *orgrermi, double *orgrerrg,
   double *orgrerpk, double *orgreltb, double *orgrelix, double *orgrelmi,
   double *orgrelrg, double *orgrelpk, double *orgimrtb, double *orgimrix,
   double *orgimrmi, double *orgimrrg, double *orgimrpk, double *orgimltb,
   double *orgimlix, double *orgimlmi, double *orgimlrg, double *orgimlpk,
   double *duprertb, double *duprerix, double *duprermi, double *duprerrg,
   double *duprerpk, double *dupreltb, double *duprelix, double *duprelmi,
   double *duprelrg, double *duprelpk, double *dupimrtb, double *dupimrix,
   double *dupimrmi, double *dupimrrg, double *dupimrpk, double *dupimltb,
   double *dupimlix, double *dupimlmi, double *dupimlrg, double *dupimlpk );
/*
 * DESCRIPTION :
 *   Makes a copy of the complex vector with real parts in ten arrays,
 *   all starting with orgre, and imaginary parts in the ten arrays
 *   prefixed with orgim, to the  complex vector with real parts in
 *   ten arrays, all starting with dupre, and imaginary parts in the ten
 *   arrays prefixed with dupim.
 *
 * REQUIRED :
 *   Space has been allocated for all vectors,
 *   to hold at least as many doubles as the value of dim.
 *
 * ON ENTRY :
 *   dim        dimension of all vectors;
 *   orgrertb   highest real parts of the original vector;
 *   orgrerix   second highest real parts of the original vector;
 *   orgrermi   third highest real parts of the original vector;
 *   orgrerrg   fourth highest real parts of the original vector;
 *   orgrerpk   fifth highest real parts of the original vector;
 *   orgreltb   fifth lowest real parts of the original vector;
 *   orgrelix   fourth lowest real parts of the original vector;
 *   orgrelmi   third lowest real parts of the original vector;
 *   orgrelrg   second lowest real parts of the original vector;
 *   orgrelpk   lowest real parts of the original vector;
 *   orgimrtb   highest imaginary parts of the original vector;
 *   orgimrix   second highest imaginary parts of the original vector;
 *   orgimrmi   third highest parts of the original vector;
 *   orgimrrg   fourth highest imaginary parts of the original vector;
 *   orgimrpk   fifth highest imaginary parts of the original vector;
 *   orgimltb   fifth lowest imaginary parts of the original vector;
 *   orgimlix   fourth lowest imaginary parts of the original vector;
 *   orgimlmi   third lowest imaginary parts of the original vector;
 *   orgimlrg   second lowest imaginary parts of the original vector;
 *   orgimlpk   lowest imaginary parts of the original vector;
 *   duprertb   space for highest real parts of a vector;
 *   duprerix   space for second highest real parts of a vector;
 *   duprermi   space for middle real parts of a vector;
 *   duprerrg   space for second lowest real parts of a vector;
 *   duprerpk   space for lowest real parts of a vector;
 *   dupreltb   space for highest real parts of a vector;
 *   duprelix   space for second highest real parts of a vector;
 *   duprelmi   space for middle real parts of a vector;
 *   duprelrg   space for second lowest real parts of a vector;
 *   duprelpk   space for lowest real parts of a vector;
 *   dupimrtb   space for highest imaginary parts of a vector;
 *   dupimrix   space for second highest imaginary parts of a vector;
 *   dupimrmi   space for third highest imaginary parts of a vector;
 *   dupimrrg   space for fourth highest imaginary parts of a vector;
 *   dupimrpk   space for fifth highest imaginary parts of a vector.
 *   dupimltb   space for fifth lowest imaginary parts of a vector;
 *   dupimlix   space for fourth lowest imaginary parts of a vector;
 *   dupimlmi   space for third lowest imaginary parts of a vector;
 *   dupimlrg   space for second lowest imaginary parts of a vector;
 *   dupimlpk   space for lowest imaginary parts of a vector.
 *
 * ON RETURN :
 *   duprertb   highest real parts of the duplicate vector;
 *   duprerix   second highest real parts of the duplicate vector;
 *   duprermi   third highest real parts of the duplicate vector;
 *   duprerrg   fourth highest real parts of the duplicate vector;
 *   duprerpk   fifth highest real parts of the duplicate vector;
 *   dupreltb   fifth lowest real parts of the duplicate vector;
 *   duprelix   fourth lowest real parts of the duplicate vector;
 *   duprelmi   third lowest real parts of the duplicate vector;
 *   duprelrg   second lowest real parts of the duplicate vector;
 *   duprelpk   lowest real parts of the duplicate vector;
 *   dupimrtb   highest imaginary parts of the duplicate vector;
 *   dupimrix   second highest imaginary parts of the duplicate vector;
 *   dupimrmi   third highest imaginary parts of the duplicate vector;
 *   dupimrrg   fourth highest imaginary parts of the duplicate vector;
 *   dupimrpk   fifth highest imaginary parts of the duplicate vector;
 *   dupimltb   fifth lowest imaginary parts of the duplicate vector;
 *   dupimlix   fourth lowest imaginary parts of the duplicate vector;
 *   dupimlmi   third lowest imaginary parts of the duplicate vector;
 *   dupimlrg   second lowest imaginary parts of the duplicate vector;
 *   dupimlpk   lowest imaginary parts of the duplicate vector. */

void CPU_norm
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk, int dim,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the complex vector given by its real parts
 *   in the first ten arrays and its imaginary parts in the next ten arrays.
 *
 * REQUIRED : All arrays have size dim.
 *
 * ON ENTRY :
 *   vrertb    highest real parts of a deca double vector;
 *   vrerix    second highest real parts of a deca double vector;
 *   vrermi    third highest real parts of a deca double vector;
 *   vrerrg    fourth highest real parts of a deca double vector;
 *   vrerpk    fifth highest real parts of a deca double vector;
 *   vreltb    fifth lowest real parts of a deca double vector;
 *   vrelix    fourth lowest real parts of a deca double vector;
 *   vrelmi    third lowest real parts of a deca double vector;
 *   vrelrg    second lowest real parts of a deca double vector;
 *   vrelpk    lowest real parts of a deca double vector;
 *   vimrtb    highest imaginary parts of a deca double vector;
 *   vimrix    second highest imaginary parts of a deca double vector;
 *   vimrmi    third highest imaginary parts of a deca double vector;
 *   vimrrg    fourth highest imaginary parts of a deca double vector;
 *   vimrpk    fifth highest imaginary parts of a deca double vector;
 *   vimltb    fifth lowest imaginary parts of a deca double vector;
 *   vimlix    fourth lowest imaginary parts of a deca double vector;
 *   vimlmi    third lowest imaginary parts of a deca double vector;
 *   vimlrg    second lowest imaginary parts of a deca double vector;
 *   vimlpk    lowest imaginary parts of a deca double vector;
 *   dim       dimension of all arrays.
 *
 * ON RETURN :
 *   normrtb   highest part of the 2-norm of the complex vector;
 *   normrix   second highest part of the 2-norm;
 *   normrmi   third highest part of the 2-norm;
 *   normrrg   fourth highest part of the 2-norm;
 *   normrpk   fifth highest part of the 2-norm;
 *   normltb   fifth lowest part of the 2-norm;
 *   normlix   fourth lowest part of the 2-norm;
 *   normlmi   third lowest part of the 2-norm;
 *   normlrg   second lowest part of the 2-norm;
 *   normlpk   lowest part of the 2-norm. */

void CPU_normalize
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk, double *vimrtb, double *vimrix,
   double *vimrmi, double *vimrrg, double *vimrpk, double *vimltb,
   double *vimlix, double *vimlmi, double *vimlrg, double *vimlpk, int dim,
   double normrtb, double normrix, double normrmi, double normrrg,
   double normrpk, double normltb, double normlix, double normlmi,
   double normlrg, double normlpk );
/*
 * DESCRIPTION :
 *   Normalizes the complex vector given by the first ten arrays for the
 *   real parts and by the next ten arrays for the imaginary parts.
 *
 * REQUIRED : All arrays have size dim.
 *
 * ON ENTRY :
 *   vrertb    highest real parts of a deca double vector;
 *   vrerix    second highest real parts of a deca double vector;
 *   vrermi    third highest real parts of a deca double vector;
 *   vrerrg    fourth highest real parts of a deca double vector;
 *   vrerpk    fifth highest real parts of a deca double vector;
 *   vreltb    fifth lowest real parts of a deca double vector;
 *   vrelix    fourth lowest real parts of a deca double vector;
 *   vrelmi    third lowest real parts of a deca double vector;
 *   vrelrg    second lowest real parts of a deca double vector;
 *   vrelpk    lowest real parts of a deca double vector;
 *   vimrtb    highest imaginary parts of a deca double vector;
 *   vimrix    second highest imaginary parts of a deca double vector;
 *   vimrmi    third highest imaginary parts of a deca double vector;
 *   vimrrg    fourth highest imaginary parts of a deca double vector;
 *   vimrpk    fifth highest imaginary parts of a deca double vector;
 *   vimltb    fifth lowest imaginary parts of a deca double vector;
 *   vimlix    fourth lowest imaginary parts of a deca double vector;
 *   vimlmi    third lowest imaginary parts of a deca double vector;
 *   vimlrg    second lowest imaginary parts of a deca double vector;
 *   vimlpk    lowest imaginary parts of a deca double vector;
 *   dim       dimension of the vector v;
 *   normrtb   highest part of the 2-norm to normalize the vector;
 *   normrix   second highest part of the 2-norm; 
 *   normrmi   third highest part of the 2-norm;
 *   normrrg   fourth highest part of the 2-norm;
 *   normrpk   fifth highest part of the 2-norm;
 *   normltb   fifth lowest part of the 2-norm; 
 *   normlix   fourth lowest part of the 2-norm;
 *   normlmi   third lowest part of the 2-norm;
 *   normlrg   second lowest part of the 2-norm;
 *   normlpk   lowest part of the 2-norm.
 *
 * ON RETURN :
 *   vrertb    highest real parts of the normalized vector;
 *   vrerix    second highest real parts of the normalized vector;
 *   vrermi    third highest real parts of the normalized vector;
 *   vrerrg    fourth highest real parts of the normalized vector;
 *   vrerpk    fifth highest real parts of the normalized vector;
 *   vreltb    fifth lowest real parts of the normalized vector;
 *   vrelix    fourth lowest real parts of the normalized vector;
 *   vrelmi    third lowest real parts of the normalized vector;
 *   vrelrg    second lowest real parts of the normalized vector;
 *   vrelpk    lowest real parts of the normalized vector;
 *   vimrtb    highest imaginary parts of the normalized vector;
 *   vimrix    second highest imaginary parts of the normalized vector;
 *   vimrmi    third highest imaginary parts of the normalized vector;
 *   vimrrg    fourth highest imaginary parts of the normalized vector;
 *   vimrpk    fifth highest imaginary parts of the normalized vector;
 *   vimltb    fifth lowest imaginary parts of the normalized vector;
 *   vimlix    fourth lowest imaginary parts of the normalized vector;
 *   vimlmi    third lowest imaginary parts of the normalized vector;
 *   vimlrg    second lowest imaginary parts of the normalized vector;
 *   vimlpk    lowest imaginary parts of the normalized vector;
 *             for the parts of the 2-norm on entry, the complex deca 
 *             double vector on return has 2-norm one. */

#endif
