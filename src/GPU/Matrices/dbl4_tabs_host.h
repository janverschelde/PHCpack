/* The file dbl4_tabs_host.h specifies functions on the host for the
 * tiled accelerated back substitution in quad double precision. */

#ifndef __dbl4_tabs_host_h__
#define __dbl4_tabs_host_h__

void CPU_dbl4_backsubs
 ( int dim, double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo );
/*
 * DESCRIPTION :
 *   Applies back substitution to solve an upper triangular system.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U
 *            and the vectors b and x;
 *   Uhihi    highest doubles of U
 *   Ulohi    second highest doubles of U
 *   Uhilo    second lowest doubles of U;
 *   Ulolo    lowest doubles of U;
 *   bhihi    highest doubles of the right hand side b;
 *   blohi    second highest doubles of the right hand side b;
 *   bhilo    second lowest doubles of the right hand side b;
 *   blolo    lowest doubles of the right hand side b;
 *   xhihi    space allocated for dim doubles;
 *   xlohi    space allocated for dim doubles;
 *   xhilo    space allocated for dim doubles;
 *   xlolo    space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhihi    highest doubles of the solution to U*x = b;
 *   xlohi    second highest doubles of the solution to U*x = b;
 *   xhilo    second lowest doubles of the solution to U*x = b;
 *   xlolo    lowest doubles of the solution to U*x = b. */

void CPU_cmplx4_backsubs
 ( int dim,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo );
/*
 * DESCRIPTION :
 *   Applies back substitution to solve an upper triangular system.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U
 *            and the vectors b and x;
 *   Urehihi  highest doubles of the real parts of U;
 *   Urelohi  second highest doubles of the real parts of U;
 *   Urehilo  second lowest doubles of the real parts of U;
 *   Urelolo  lowest doubles of the real parts of U;
 *   Uimhihi  highest doubles of the imaginary parts of U;
 *   Uimlohi  second highest doubles of the imaginary parts of U;
 *   Uimhilo  second lowest doubles of the imaginary parts of U;
 *   Uimlolo  lowest doubles of the imaginary parts of U;
 *   brehihi  highest doubles of the real parts of b;
 *   brelohi  second highest doubles of the real parts of b;
 *   brehilo  second lowest doubles of the real parts of b;
 *   brelolo  lowest doubles of the real parts of b;
 *   bimhihi  highest doubles of the imaginary parts of b;
 *   bimlohi  second highest doubles of the imaginary parts of b;
 *   bimhilo  second lowest doubles of the imaginary parts of b;
 *   bimlolo  lowest doubles of the imaginary parts of b;
 *   xrehihi  space allocated for dim doubles;
 *   xrelolo  space allocated for dim doubles;
 *   ximhihi  space allocated for dim doubles;
 *   ximlolo  space allocated for dim doubles.
 *
 * ON RETURN :
 *   xrehihi  highest doubles of the real parts of the solution x;
 *   xrelohi  second highest doubles of the real parts of the solution x;
 *   xrehilo  second lowest doubles of the real parts of the solution x;
 *   xrelolo  lowest doubles of the real parts of the solution x;
 *   ximhihi  highest doubles of the imaginary parts of the solution x;
 *   ximlohi  second highest doubles of the imaginary parts of the solution x;
 *   ximhilo  second lowest doubles of the imaginary parts of the solution x;
 *   ximlolo  lowest doubles of the imaginary parts of the solution x. */

void CPU_dbl4_upper_inverse
 ( int dim,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double **invUhihi, double **invUlohi, double **invUhilo, double **invUlolo,
   double *lapsec );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhihi    highest doubles of U;
 *   Ulohi    second highest doubles of U;
 *   Uhilo    second lowest doubles of U;
 *   Ulolo    lowest doubles of U;
 *   invUhihi has space allocated for a matrix of dimension dim;
 *   invUlohi has space allocated for a matrix of dimension dim;
 *   invUhilo has space allocated for a matrix of dimension dim;
 *   invUlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhihi has the highest doubles of the inverse of the matrix U;
 *   invUlohi has the second highest doubles of the inverse of the matrix U;
 *   invUhilo has the second lowest doubles of the inverse of the matrix U;
 *   invUlolo has the lowest doubles of the inverse of the matrix U;
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx4_upper_inverse
 ( int dim,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double **invUrehihi, double **invUrelohi,
   double **invUrehilo, double **invUrelolo,
   double **invUimhihi, double **invUimlohi,
   double **invUimhilo, double **invUimlolo, double *lapsec );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Urehihi  highest doubles of the real parts of U;
 *   Urelohi  second highest doubles of the real parts of U;
 *   Urehilo  second lowest doubles of the real parts of U;
 *   Urelolo  lowest doubles of the real parts of U;
 *   Uimhihi  second highest doubles of the imaginary parts of U;
 *   Uimlohi  second highest doubles of the imaginary parts of U;
 *   Uimhilo  lowest doubles of the imaginary parts of U;
 *   Uimlolo  lowest doubles of the imaginary parts of U;
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
 *   invUrehihi has the highest doubles of the real parts of the inverse;
 *   invUrelohi has the second highest doubles of the real parts
 *            of the inverse;
 *   invUrehilo has the second lowest doubles of the real parts
 *            of the inverse;
 *   invUrehilo has the lowest doubles of the real parts of the inverse;
 *   invUimhihi has the highest doubles of the imaginary parts of the inverse;
 *   invUimlohi has the second highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimhilo has the second lowest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlolo has the lowest doubles of the imaginary parts of the inverse;
 *   lapsec   elapsed time in seconds. */

void CPU_dbl4_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *lapsec );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Uhihi    highest doubles of U;
 *   Ulohi    second highest doubles of U;
 *   Uhilo    second lowest doubles of U;
 *   Ulolo    lowest doubles of U;
 *   bhihi    highest doubles of b;
 *   blohi    second highest doubles of b;
 *   bhilo    second lowest doubles of b;
 *   blolo    lowest doubles of b;
 *   xhihi    space allocated for dim doubles;
 *   xlohi    space allocated for dim doubles;
 *   xhilo    space allocated for dim doubles;
 *   xlolo    space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhihi    highest doubles of the solution to U*x = b;
 *   xlohi    second highest doubles of the solution to U*x = b;
 *   xhilo    second lowest doubles of the solution to U*x = b;
 *   xlolo    lowest doubles of the solution to U*x = b;
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx4_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *lapsec );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Urehihi  highest doubles of the real parts of U;
 *   Urelohi  second highest doubles of the real parts of U;
 *   Urehilo  second lowest doubles of the real parts of U;
 *   Urelolo  lowest doubles of the real parts of U;
 *   Uimhihi  highest doubles of the imaginary parts of U;
 *   Uimlohi  second highest doubles of the imaginary parts of U;
 *   Uimhilo  second lowest doubles of the imaginary parts of U;
 *   Uimlolo  lowest doubles of the imaginary parts of U;
 *   brehihi  highest doubles of the real parts of b;
 *   brelohi  second highest doubles of the real parts of b;
 *   brehilo  second lowest doubles of the real parts of b;
 *   brelolo  lowest doubles of the real parts of b;
 *   bimhihi  highest doubles of the imaginary parts of b;
 *   bimlohi  second highest doubles of the imaginary parts of b;
 *   bimhilo  second lowest doubles of the imaginary parts of b;
 *   bimlolo  lowest doubles of the imaginary parts of b;
 *   xrehihi  space allocated for dim doubles;
 *   xrelohi  space allocated for dim doubles;
 *   xrehilo  space allocated for dim doubles;
 *   xrelolo  space allocated for dim doubles;
 *   ximhihi  space allocated for dim doubles;
 *   ximlohi  space allocated for dim doubles;
 *   ximhilo  space allocated for dim doubles;
 *   ximlolo  space allocated for dim doubles.
 *
 * ON RETURN :
 *   Urehihi  the diagonal tiles contain the highest doubles of
 *            the real parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Urelohi  the diagonal tiles contain the second highest doubles of
 *            the real parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Urehilo  the diagonal tiles contain the second lowest doubles of
 *            the real parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Urelolo  the diagonal tiles contain the lowest doubles of
 *            the real parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Uimhihi  the diagonal tiles contain the highest doubles of
 *            the imaginary parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Uimlohi  the diagonal tiles contain the second highest doubles of
 *            the imaginary parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Uimhilo  the diagonal tiles contain the second lowest doubles of
 *            the imaginary parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Uimlolo  the diagonal tiles contain the lowest doubles of
 *            the imaginary parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   xrehi    highest doubles of the real parts of the solution x;
 *   xrehi    second highest doubles of the real parts of the solution x;
 *   xrelo    second lowest doubles of the real parts of the solution x;
 *   xrelo    lowest doubles of the real parts of the solution x;
 *   ximhi    highest doubles of the imaginary parts of the solution x;
 *   ximhi    second highest doubles of the imaginary parts of the solution x;
 *   ximlo    second lowest doubles of the imaginary parts of the solution x;
 *   ximlo    lowest doubles of the imaginary parts of the solution x;
 *   lapsec   elapsed time in seconds. */

#endif
