// The file cmplx8_norm_host.h specifies functions to compute
// the 2-norm of a complex vector and to normalize a complex vector,
// in octo double precision.

#ifndef __cmplx8_norm_host_h__
#define __cmplx8_norm_host_h__

void make_copy
 ( int dim,
   double *orgrehihihi, double *orgrelohihi, double *orgrehilohi,
   double *orgrelolohi, double *orgrehihilo, double *orgrelohilo,
   double *orgrehilolo, double *orgrelololo,
   double *orgimhihihi, double *orgimlohihi, double *orgimhilohi,
   double *orgimlolohi, double *orgimhihilo, double *orgimlohilo,
   double *orgimhilolo, double *orgimlololo,
   double *duprehihihi, double *duprelohihi, double *duprehilohi,
   double *duprelolohi, double *duprehihilo, double *duprelohilo,
   double *duprehilolo, double *duprelololo,
   double *dupimhihihi, double *dupimlohihi, double *dupimhilohi,
   double *dupimlolohi, double *dupimhihilo, double *dupimlohilo,
   double *dupimhilolo, double *dupimlololo );
/*
 * DESCRIPTION :
 *   Makes a copy of the complex vector org with real parts in eight arrays
 *   and imaginary parts in eight arrays to the complex vector dup with
 *   corresponding two sequences of eight arrays.
 *
 * REQUIRED :
 *   Space has been allocated for all vectors,
 *   to hold at least as many doubles as the value of dim.
 *
 * ON ENTRY :
 *   dim           dimension of all vectors;
 *   orgrehihihi   highest real parts of the original vector;
 *   orgrelohihi   second highest real parts of the original vector;
 *   orgrehilohi   third highest real parts of the original vector;
 *   orgrelolohi   fourth highest real parts of the original vector;
 *   orgrehihilo   fourth lowest real parts of the original vector;
 *   orgrelohilo   third lowest real parts of the original vector;
 *   orgrehilolo   second lowest real parts of the original vector;
 *   orgrelololo   lowest real parts of the original vector;
 *   orgimhihihi   highest imaginary parts of the original vector;
 *   orgimlohihi   second highest imaginary parts of the original vector;
 *   orgimhilohi   third highest imaginary parts of the original vector;
 *   orgimlolohi   fourth highest imaginary parts of the original vector;
 *   orgimhihilo   fourth lowest imaginary parts of the original vector;
 *   orgimlohilo   third lowest imaginary parts of the original vector;
 *   orgimhilolo   second lowest imaginary parts of the original vector;
 *   orgimlololo   lowest imaginary parts of the original vector;
 *   duprehihihi   space for highest real parts of a vector;
 *   duprelohihi   space for second highest real parts of a vector;
 *   duprehilohi   space for third highest real parts of a vector;
 *   duprelolohi   space for fourth highest real parts of a vector;
 *   duprehihilo   space for fourth lowest real parts of a vector;
 *   duprelohilo   space for third lowest real parts of a vector;
 *   duprehilolo   space for second lowest real parts of a vector;
 *   duprelololo   space for lowest real parts of a vector;
 *   dupimhihihi   space for highest imaginary parts of a vector;
 *   dupimlohihi   space for second highest imaginary parts of a vector;
 *   dupimhilohi   space for third highest imaginary parts of a vector;
 *   dupimlolohi   space for fourth highest imaginary parts of a vector.
 *   dupimhihilo   space for fourth lowest imaginary parts of a vector;
 *   dupimlohilo   space for third lowest imaginary parts of a vector;
 *   dupimhilolo   space for second lowest imaginary parts of a vector;
 *   dupimlololo   space for lowest imaginary parts of a vector.
 *
 * ON RETURN :
 *   duprehihihi   highest real parts of the duplicate vector;
 *   duprelohihi   second highest real parts of the duplicate vector;
 *   duprehilohi   third highest real parts of the duplicate vector;
 *   duprelolohi   fourth highest real parts of the duplicate vector;
 *   duprehihilo   fourth lowest real parts of the duplicate vector;
 *   duprelohilo   third lowest real parts of the duplicate vector;
 *   duprehilolo   second lowest real parts of the duplicate vector;
 *   duprelololo   lowest real parts of the duplicate vector;
 *   dupimhihihi   highest imaginary parts of the duplicate vector;
 *   dupimlohihi   second highest imaginary parts of the duplicate vector;
 *   dupimhilohi   third highest imaginary parts of the duplicate vector;
 *   dupimlolohi   fourth highest imaginary parts of the duplicate vector;
 *   dupimhihilo   fourth lowest imaginary parts of the duplicate vector;
 *   dupimlohilo   third lowest imaginary parts of the duplicate vector;
 *   dupimhilolo   second lowest imaginary parts of the duplicate vector;
 *   dupimlololo   lowest imaginary parts of the duplicate vector. */

void CPU_norm
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi,
   double *vrelolohi, double *vrehihilo, double *vrelohilo,
   double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi,
   double *vimlolohi, double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo,
   int dim,
   double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the complex vector given by its real parts
 *   in eight arrays and its imaginary parts in eight arrays.
 *
 * REQUIRED : All arrays have size dim.
 *
 * ON ENTRY :
 *   vrehihihi   highest real parts of a octo double vector;
 *   vrelohihi   second highest real parts of a octo double vector;
 *   vrehilohi   third highest real parts of a octo double vector;
 *   vrelolohi   fourth highest real parts of a octo double vector;
 *   vrehihilo   fourth lowest real parts of a octo double vector;
 *   vrelohilo   third lowest real parts of a octo double vector;
 *   vrehilolo   second lowest real parts of a octo double vector;
 *   vrelololo   lowest real parts of a octo double vector;
 *   vimhihihi   highest imaginary parts of a octo double vector;
 *   vimlohihi   second highest imaginary parts of a octo double vector;
 *   vimhilohi   third highest imaginary parts of a octo double vector;
 *   vimlolohi   fourth highest imaginary parts of a octo double vector;
 *   vimhihilo   fourth lowest imaginary parts of a octo double vector;
 *   vimlohilo   third lowest imaginary parts of a octo double vector;
 *   vimhilolo   second lowest imaginary parts of a octo double vector;
 *   vimlololo   lowest imaginary parts of a octo double vector;
 *   dim         dimension of all arrays.
 *
 * ON RETURN :
 *   normhihihi  highest part of the 2-norm of the complex vector;
 *   normlohihi  second highest part of the 2-norm of the complex vector;
 *   normhilohi  third highest part of the 2-norm of the complex vector;
 *   normlolohi  fourth highest part of the 2-norm of the complex vector;
 *   normhihilo  fourth lowest part of the 2-norm of the complex vector;
 *   normlohilo  third lowest part of the 2-norm of the complex vector;
 *   normhilolo  second lowest part of the 2-norm of the complex vector;
 *   normlololo  lowest part of the 2-norm of the complex vector. */

void CPU_normalize
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi,
   double *vrelolohi, double *vrehihilo, double *vrelohilo,
   double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi,
   double *vimlolohi, double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo,
   int dim,
   double normhihihi, double normlohihi, double normhilohi,
   double normlolohi, double normhihilo, double normlohilo,
   double normhilolo, double normlololo );
/*
 * DESCRIPTION :
 *   Normalizes the complex vector given by eight arrays for the
 *   real parts and by eight arrays for the imaginary parts.
 *
 * REQUIRED : All arrays have size dim.
 *
 * ON ENTRY :
 *   vrehihihi   highest real parts of a octo double vector;
 *   vrelohihi   second highest real parts of a octo double vector;
 *   vrehilohi   third highest real parts of a octo double vector;
 *   vrelolohi   fourth highest real parts of a octo double vector;
 *   vrehihilo   fourth lowest real parts of a octo double vector;
 *   vrelohilo   third lowest real parts of a octo double vector;
 *   vrehilolo   second lowest real parts of a octo double vector;
 *   vrelololo   lowest real parts of a octo double vector;
 *   vimhihihi   highest imaginary parts of a octo double vector;
 *   vimlohihi   second highest imaginary parts of a octo double vector;
 *   vimhilohi   third highest imaginary parts of a octo double vector;
 *   vimlolohi   fourth highest imaginary parts of a octo double vector;
 *   vimhihilo   fourth lowest imaginary parts of a octo double vector;
 *   vimlohilo   third lowest imaginary parts of a octo double vector;
 *   vimhilolo   second lowest imaginary parts of a octo double vector;
 *   vimlololo   lowest imaginary parts of a octo double vector;
 *   dim         dimension of the vector v;
 *   normhihihi  highest part of the norm to normalize the vector;
 *   normlohihi  second highest part of the norm to normalize the vector;
 *   normhilohi  third highest part of the norm to normalize the vector;
 *   normlolohi  fourth highest part of the norm to normalize the vector.
 *   normhihilo  fourth lowest part of the norm to normalize the vector;
 *   normlohilo  third lowest part of the norm to normalize the vector;
 *   normhilolo  second lowest part of the norm to normalize the vector;
 *   normlololo  lowest part of the norm to normalize the vector.
 *
 * ON RETURN :
 *   vrehihihi   highest real parts of the normalized vector;
 *   vrelohihi   second highest real parts of the normalized vector;
 *   vrehilohi   third highest real parts of the normalized vector;
 *   vrelolohi   fourth highest real parts of the normalized vector;
 *   vrehihilo   fourth lowest real parts of the normalized vector;
 *   vrelohilo   third lowest real parts of the normalized vector;
 *   vrehilolo   second lowest real parts of the normalized vector;
 *   vrelololo   lowest real parts of the normalized vector;
 *   vimhihihi   highest imaginary parts of the normalized vector;
 *   vimlohihi   second highest imaginary parts of the normalized vector;
 *   vimhilohi   third highest imaginary parts of the normalized vector;
 *   vimlolohi   fourth highest imaginary parts of the normalized vector;
 *   vimhihilo   fourth lowest imaginary parts of the normalized vector;
 *   vimlohilo   third lowest imaginary parts of the normalized vector;
 *   vimhilolo   second lowest imaginary parts of the normalized vector;
 *   vimlololo   lowest imaginary parts of the normalized vector;
 *               if on entry, the norm* numbers equal the 2-norm,
 *               then the complex octo double vector on return 
 *               has 2-norm one. */

#endif
