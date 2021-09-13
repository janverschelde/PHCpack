/* The file dbl8_tabs_host.h specifies functions on the host for the
 * tiled accelerated back substitution in octo double precision. */

#ifndef __dbl8_tabs_host_h__
#define __dbl8_tabs_host_h__

void CPU_dbl8_backsubs
 ( int dim,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo );
/*
 * DESCRIPTION :
 *   Applies back substitution to solve an upper triangular system.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U
 *            and the vectors b and x;
 *   Uhihihi  are the highest doubles of U
 *   Ulohihi  are the second highest doubles of U
 *   Uhilohi  are the third highest doubles of U
 *   Ulolohi  are the fourth highest doubles of U
 *   Uhihilo  are the fourth lowest doubles of U;
 *   Ulohilo  are the third lowest doubles of U;
 *   Uhilolo  are the second lowest doubles of U;
 *   Ulololo  are the lowest doubles of U;
 *   bhihihi  are the highest doubles of b;
 *   blohihi  are the second highest doubles of b;
 *   bhilohi  are the third highest doubles of b;
 *   blolohi  are the fourth highest doubles of b;
 *   bhihilo  are the fourth lowest doubles of b;
 *   blohilo  are the third lowest doubles of b;
 *   bhilolo  are the second lowest doubles of b;
 *   blololo  are the lowest doubles of b;
 *   xhihihi  has space allocated for dim doubles;
 *   xlohihi  has space allocated for dim doubles;
 *   xhilohi  has space allocated for dim doubles;
 *   xlolohi  has space allocated for dim doubles;
 *   xhihilo  has space allocated for dim doubles;
 *   xlohilo  has space allocated for dim doubles;
 *   xhilolo  has space allocated for dim doubles;
 *   xlololo  has space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhihihi  are the highest doubles of the solution x;
 *   xlohihi  are the second highest doubles of x;
 *   xhilohi  are the third highest doubles of x;
 *   xlolohi  are the fourth highest doubles of x;
 *   xhihilo  are the fourth lowest doubles of x;
 *   xlohilo  are the third lowest doubles of x;
 *   xhilolo  are the second lowest doubles of x;
 *   xlololo  are the lowest doubles of x. */

void CPU_cmplx8_backsubs
 ( int dim,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo );
/*
 * DESCRIPTION :
 *   Applies back substitution to solve an upper triangular system.
 *
 * ON ENTRY :
 *   dim       dimension of the upper triangular matrix U
 *             and the vectors b and x;
 *   Urehihihi are the highest doubles of the real parts of U;
 *   Urelohihi are the second highest doubles of the real parts of U;
 *   Urehilohi are the third highest doubles of the real parts of U;
 *   Urelolohi are the fourth highest doubles of the real parts of U;
 *   Urehihilo are the fourth lowest doubles of the real parts of U;
 *   Urelohilo are the third lowest doubles of the real parts of U;
 *   Urehilolo are the second lowest doubles of the real parts of U;
 *   Urelololo are the lowest doubles of the real parts of U;
 *   Uimhihihi are the highest doubles of the imaginary parts of U;
 *   Uimlohihi are the second highest doubles of the imaginary parts of U;
 *   Uimhilohi are the third highest doubles of the imaginary parts of U;
 *   Uimlolohi are the fourth highest doubles of the imaginary parts of U;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts of U;
 *   Uimlohilo are the third lowest doubles of the imaginary parts of U;
 *   Uimhilolo are the second lowest doubles of the imaginary parts of U;
 *   Uimlololo are the lowest doubles of the imaginary parts of U;
 *   brehihihi are the highest doubles of the real parts of b;
 *   brelohihi are the second highest doubles of the real parts of b;
 *   brehilohi are the third highest doubles of the real parts of b;
 *   brelolohi are the fourth highest doubles of the real parts of b;
 *   brehihilo are the fourth lowest doubles of the real parts of b;
 *   brelohilo are the third lowest doubles of the real parts of b;
 *   brehilolo are the second lowest doubles of the real parts of b;
 *   brelololo are the lowest doubles of the real parts of b;
 *   bimhihihi are the highest doubles of the imaginary parts of b;
 *   bimlohihi are the second highest doubles of the imaginary parts of b;
 *   bimhilohi are the third highest doubles of the imaginary parts of b;
 *   bimlolohi are the fourth highest doubles of the imaginary parts of b;
 *   bimhihilo are the fourth lowest doubles of the imaginary parts of b;
 *   bimlohilo are the third lowest doubles of the imaginary parts of b;
 *   bimhilolo are the second lowest doubles of the imaginary parts of b;
 *   bimlololo are the lowest doubles of the imaginary parts of b;
 *   xrehihihi has space allocated for dim doubles;
 *   xrelohihi has space allocated for dim doubles;
 *   xrehilohi has space allocated for dim doubles;
 *   xrelolohi has space allocated for dim doubles;
 *   xrehihilo has space allocated for dim doubles;
 *   xrelohilo has space allocated for dim doubles;
 *   xrehilolo has space allocated for dim doubles;
 *   xrelololo has space allocated for dim doubles;
 *   ximhihihi has space allocated for dim doubles;
 *   ximlohihi has space allocated for dim doubles;
 *   ximhilohi has space allocated for dim doubles;
 *   ximlolohi has space allocated for dim doubles.
 *   ximhihilo has space allocated for dim doubles;
 *   ximlohilo has space allocated for dim doubles.
 *   ximhilolo has space allocated for dim doubles;
 *   ximlololo has space allocated for dim doubles.
 *
 * ON RETURN :
 *   xrehihihi are the highest doubles of the real parts of the solution x;
 *   xrelohihi are the second highest doubles of the real parts of x;
 *   xrehilohi are the third highest doubles of the real parts of x;
 *   xrelolohi are the fourth highest doubles of the real parts of x;
 *   xrehihilo are the fourth lowest doubles of the real parts of x;
 *   xrelohilo are the third lowest doubles of the real parts of x;
 *   xrehilolo are the second lowest doubles of the real parts of x;
 *   xrelololo are the lowest doubles of the real parts of x;
 *   ximhihihi are the highest doubles of the imaginary parts of x;
 *   ximlohihi are the second highest doubles of the imaginary parts of x;
 *   ximhilohi are the third highest doubles of the imaginary parts of x;
 *   ximlolohi are the fourth highest doubles of the imaginary parts of x;
 *   ximhihilo are the fourth lowest doubles of the imaginary parts of x;
 *   ximlohilo are the third lowest doubles of the imaginary parts of x;
 *   ximhilolo are the second lowest doubles of the imaginary parts of x;
 *   ximlololo are the lowest doubles of the imaginary parts of x. */

void CPU_dbl8_upper_inverse
 ( int dim,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double **invUhihihi, double **invUlohihi,
   double **invUhilohi, double **invUlolohi,
   double **invUhihilo, double **invUlohilo,
   double **invUhilolo, double **invUlololo, double *lapsec );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhihihi  are the highest doubles of U
 *   Ulohihi  are the second highest doubles of U
 *   Uhilohi  are the third highest doubles of U
 *   Ulolohi  are the fourth highest doubles of U
 *   Uhihilo  are the fourth lowest doubles of U;
 *   Ulohilo  are the third lowest doubles of U;
 *   Uhilolo  are the second lowest doubles of U;
 *   Ulololo  are the lowest doubles of U;
 *   invUhihihi has space allocated for a matrix of dimension dim;
 *   invUlohihi has space allocated for a matrix of dimension dim;
 *   invUhilolo has space allocated for a matrix of dimension dim;
 *   invUlololo has space allocated for a matrix of dimension dim;
 *   invUhihihi has space allocated for a matrix of dimension dim;
 *   invUlohihi has space allocated for a matrix of dimension dim;
 *   invUhilolo has space allocated for a matrix of dimension dim;
 *   invUlololo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhihihi has the highest doubles of the inverse of U;
 *   invUlohihi has the second highest doubles of the inverse of U;
 *   invUhilohi has the third highest doubles of the inverse of U;
 *   invUlolohi has the fourth highest doubles of the inverse of U;
 *   invUhihilo has the fourth lowest doubles of the inverse of U;
 *   invUlohilo has the third lowest doubles of the inverse of U;
 *   invUhilolo has the second lowest doubles of the inverse of U;
 *   invUlololo has the lowest doubles of the inverse of U;
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx8_upper_inverse
 ( int dim,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double **invUrehihihi, double **invUrelohihi,
   double **invUrehilohi, double **invUrelolohi,
   double **invUrehihilo, double **invUrelohilo,
   double **invUrehilolo, double **invUrelololo,
   double **invUimhihihi, double **invUimlohihi,
   double **invUimhilohi, double **invUimlolohi,
   double **invUimhihilo, double **invUimlohilo,
   double **invUimhilolo, double **invUimlololo, double *lapsec );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *
 * ON ENTRY :
 *   dim       dimension of the upper triangular matrix U;
 *   Urehihihi are the highest doubles of the real parts of U;
 *   Urelohihi are the second highest doubles of the real parts of U;
 *   Urehilohi are the third highest doubles of the real parts of U;
 *   Urelolohi are the fourth highest doubles of the real parts of U;
 *   Urehihilo are the fourth lowest doubles of the real parts of U;
 *   Urelohilo are the third lowest doubles of the real parts of U;
 *   Urehilolo are the second lowest doubles of the real parts of U;
 *   Urelololo are the lowest doubles of the real parts of U;
 *   Uimhihihi are the highest doubles of the imaginary parts of U;
 *   Uimlohihi are the second highest doubles of the imaginary parts of U;
 *   Uimhilohi are the third highest doubles of the imaginary parts of U;
 *   Uimlolohi are the fourth highest doubles of the imaginary parts of U;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts of U;
 *   Uimlohilo are the third lowest doubles of the imaginary parts of U;
 *   Uimhilolo are the second lowest doubles of the imaginary parts of U;
 *   Uimlololo are the lowest doubles of the imaginary parts of U;
 *   invUrehihi has space allocated for a matrix of dimension dim;
 *   invUrelohi has space allocated for a matrix of dimension dim;
 *   invUrehilo has space allocated for a matrix of dimension dim;
 *   invUrelolo has space allocated for a matrix of dimension dim;
 *   invUimhihi has space allocated for a matrix of dimension dim;
 *   invUimlohi has space allocated for a matrix of dimension dim;
 *   invUimhilo has space allocated for a matrix of dimension dim;
 *   invUimlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehihihi has the highest doubles of the real parts of the inverse;
 *   invUrelohihi has the second highest doubles of the real parts
 *            of the inverse;
 *   invUrehilohi has the third highest doubles of the real parts
 *            of the inverse;
 *   invUrelolohi has the fourth highest doubles of the real parts
 *            of the inverse;
 *   invUrehihilo has the fourth lowest doubles of the real parts
 *            of the inverse;
 *   invUrelohilo has the third lowest doubles of the real parts
 *            of the inverse;
 *   invUrehilolo has the second lowest doubles of the real parts
 *            of the inverse;
 *   invUrelololo has the lowest doubles of the real parts of the inverse;
 *   invUimhihihi has the highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlohihi has the second highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimhilohi has the third highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlolohi has the fourth highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimhihilo has the fourth lowest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlohilo has the third lowest doubles of the imaginary parts
 *            of the inverse;
 *   invUimhilolo has the second lowest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlololo has the lowest doubles of the imaginary parts
 *            of the inverse;
 *   lapsec   elapsed time in seconds. */

void CPU_dbl8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *lapsec );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Uhihihi  are the highest doubles of U
 *   Ulohihi  are the second highest doubles of U
 *   Uhilohi  are the third highest doubles of U
 *   Ulolohi  are the fourth highest doubles of U
 *   Uhihilo  are the fourth lowest doubles of U;
 *   Ulohilo  are the third lowest doubles of U;
 *   Uhilolo  are the second lowest doubles of U;
 *   Ulololo  are the lowest doubles of U;
 *   bhihihi  are the highest doubles of b;
 *   blohihi  are the second highest doubles of b;
 *   bhilohi  are the third highest doubles of b;
 *   blolohi  are the fourth highest doubles of b;
 *   bhihilo  are the fourth lowest doubles of b;
 *   blohilo  are the third lowest doubles of b;
 *   bhilolo  are the second lowest doubles of b;
 *   blololo  are the lowest doubles of b;
 *   xhihihi  has space allocated for dim doubles;
 *   xlohihi  has space allocated for dim doubles;
 *   xhilohi  has space allocated for dim doubles;
 *   xlolohi  has space allocated for dim doubles;
 *   xhihilo  has space allocated for dim doubles;
 *   xlohilo  has space allocated for dim doubles;
 *   xhilolo  has space allocated for dim doubles;
 *   xlololo  has space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhihihi  are the highest doubles of the solution x;
 *   xlohihi  are the second highest doubles of x;
 *   xhilohi  are the third highest doubles of x;
 *   xlolohi  are the fourth highest doubles of x;
 *   xhihilo  are the fourth lowest doubles of x;
 *   xlohilo  are the third lowest doubles of x;
 *   xhilolo  are the second lowest doubles of x;
 *   xlololo  are the lowest doubles of x;
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *lapsec );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim       dimension of the upper triangular matrix U;
 *   szt       size of each tile;
 *   nbt       number of tiles, dim = szt*nbt;
 *   Urehihihi are the highest doubles of the real parts of U;
 *   Urelohihi are the second highest doubles of the real parts of U;
 *   Urehilohi are the third highest doubles of the real parts of U;
 *   Urelolohi are the fourth highest doubles of the real parts of U;
 *   Urehihilo are the fourth lowest doubles of the real parts of U;
 *   Urelohilo are the third lowest doubles of the real parts of U;
 *   Urehilolo are the second lowest doubles of the real parts of U;
 *   Urelololo are the lowest doubles of the real parts of U;
 *   Uimhihihi are the highest doubles of the imaginary parts of U;
 *   Uimlohihi are the second highest doubles of the imaginary parts of U;
 *   Uimhilohi are the third highest doubles of the imaginary parts of U;
 *   Uimlolohi are the fourth highest doubles of the imaginary parts of U;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts of U;
 *   Uimlohilo are the third lowest doubles of the imaginary parts of U;
 *   Uimhilolo are the second lowest doubles of the imaginary parts of U;
 *   Uimlololo are the lowest doubles of the imaginary parts of U;
 *   brehihihi are the highest doubles of the real parts of b;
 *   brelohihi are the second highest doubles of the real parts of b;
 *   brehilohi are the third highest doubles of the real parts of b;
 *   brelolohi are the fourth highest doubles of the real parts of b;
 *   brehihilo are the fourth lowest doubles of the real parts of b;
 *   brelohilo are the third lowest doubles of the real parts of b;
 *   brehilolo are the second lowest doubles of the real parts of b;
 *   brelololo are the lowest doubles of the real parts of b;
 *   bimhihihi are the highest doubles of the imaginary parts of b;
 *   bimlohihi are the second highest doubles of the imaginary parts of b;
 *   bimhilohi are the third highest doubles of the imaginary parts of b;
 *   bimlolohi are the fourth highest doubles of the imaginary parts of b;
 *   bimhihilo are the fourth lowest doubles of the imaginary parts of b;
 *   bimlohilo are the third lowest doubles of the imaginary parts of b;
 *   bimhilolo are the second lowest doubles of the imaginary parts of b;
 *   bimlololo are the lowest doubles of the imaginary parts of b;
 *   xrehihihi has space allocated for dim doubles;
 *   xrelohihi has space allocated for dim doubles;
 *   xrehilohi has space allocated for dim doubles;
 *   xrelolohi has space allocated for dim doubles;
 *   xrehihilo has space allocated for dim doubles;
 *   xrelohilo has space allocated for dim doubles;
 *   xrehilolo has space allocated for dim doubles;
 *   xrelololo has space allocated for dim doubles;
 *   ximhihihi has space allocated for dim doubles;
 *   ximlohihi has space allocated for dim doubles;
 *   ximhilohi has space allocated for dim doubles;
 *   ximlolohi has space allocated for dim doubles.
 *   ximhihilo has space allocated for dim doubles;
 *   ximlohilo has space allocated for dim doubles.
 *   ximhilolo has space allocated for dim doubles;
 *   ximlololo has space allocated for dim doubles.
 *
 * ON RETURN :
 *   xrehihihi are the highest doubles of the real parts of the solution x;
 *   xrelohihi are the second highest doubles of the real parts of x;
 *   xrehilohi are the third highest doubles of the real parts of x;
 *   xrelolohi are the fourth highest doubles of the real parts of x;
 *   xrehihilo are the fourth lowest doubles of the real parts of x;
 *   xrelohilo are the third lowest doubles of the real parts of x;
 *   xrehilolo are the second lowest doubles of the real parts of x;
 *   xrelololo are the lowest doubles of the real parts of x;
 *   ximhihihi are the highest doubles of the imaginary parts of x;
 *   ximlohihi are the second highest doubles of the imaginary parts of x;
 *   ximhilohi are the third highest doubles of the imaginary parts of x;
 *   ximlolohi are the fourth highest doubles of the imaginary parts of x;
 *   ximhihilo are the fourth lowest doubles of the imaginary parts of x;
 *   ximlohilo are the third lowest doubles of the imaginary parts of x;
 *   ximhilolo are the second lowest doubles of the imaginary parts of x;
 *   ximlololo are the lowest doubles of the imaginary parts of x;
 *   lapsec   elapsed time in seconds. */

#endif
