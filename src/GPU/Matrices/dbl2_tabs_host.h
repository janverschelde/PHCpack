/* The file dbl2_tabs_host.h specifies functions on the host for the
 * tiled accelerated back substitution in double double precision. */

#ifndef __dbl2_tabs_host_h__
#define __dbl2_tabs_host_h__

void CPU_dbl2_backsubs
 ( int dim, double **Uhi, double **Ulo, double *bhi, double *blo,
   double *xhi, double *xlo );
/*
 * DESCRIPTION :
 *   Applies back substitution to solve an upper triangular system.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U
 *            and the vectors b and x;
 *   Uhi      high doubles of an upper triangular matrix of dimension dim;
 *   Ulo      low doubles of an upper triangular matrix of dimension dim;
 *   bhi      high doubles of the right hand side of the linear system;
 *   blo      low doubles of the right hand side of the linear system;
 *   xhi      space allocated for dim doubles;
 *   xlo      space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhi      high doubles of the solution to the system U*x = b;
 *   xlo      low doubles of the solution to the system U*x = b. */

void CPU_cmplx2_backsubs
 ( int dim, double **Urehi, double **Urelo,
            double **Uimhi, double **Uimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo );
/*
 * DESCRIPTION :
 *   Applies back substitution to solve an upper triangular system.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U
 *            and the vectors b and x;
 *   Urehi    high doubles of the real parts of U;
 *   Urelo    low doubles of the real parts of U;
 *   Uimhi    high doubles of the imaginary parts of U;
 *   Uimlo    low doubles of the imaginary parts of U;
 *   brehi    high doubles of the real parts of b;
 *   brelo    low doubles of the real parts of b;
 *   bimhi    high doubles of the imaginary parts of b;
 *   bimlo    low doubles of the imaginary parts of b;
 *   xrehi    space allocated for dim doubles;
 *   xrelo    space allocated for dim doubles;
 *   ximhi    space allocated for dim doubles;
 *   ximlo    space allocated for dim doubles.
 *
 * ON RETURN :
 *   xrehi    high doubles of the real parts of the solution x;
 *   xrelo    low doubles of the real parts of the solution x;
 *   ximhi    high doubles of the imaginary parts of the solution x;
 *   ximlo    low doubles of the imaginary parts of the solution x. */

void CPU_dbl2_upper_inverse
 ( int dim, double **Uhi, double **Ulo, double **invUhi, double **invUlo,
   double *lapsec );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhi      high doubles of an upper triangular matrix of dimension dim;
 *   Ulo      low doubles of an upper triangular matrix of dimension dim;
 *   invUhi   space allocated for a matrix of dimension dim;
 *   invUlo   space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhi   high doubles of the inverse of the matrix U;
 *   invUlo   low doubles of the inverse of the matrix U;
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx2_upper_inverse
 ( int dim, double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double **invUrehi, double **invUrelo,
   double **invUimhi, double **invUimlo, double *lapsec );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Urehi    high doubles of the real parts of U;
 *   Urelo    low doubles of the real parts of U;
 *   Uimhi    high doubles of the imaginary parts of U;
 *   Uimlo    low doubles of the imaginary parts of U;
 *   invUrehi has space allocated for a matrix of dimension dim;
 *   invUrelo has space allocated for a matrix of dimension dim;
 *   invUimhi has space allocated for a matrix of dimension dim;
 *   invUimlo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehi has the high doubles of the real parts of the inverse;
 *   invUrelo has the low doubles of the real parts of the inverse;
 *   invUimhi has the high doubles of the imaginary parts of the inverse;
 *   invUimlo has the low doubles of the imaginary parts of the inverse;
 *   lapsec   elapsed time in seconds. */

void CPU_dbl2_matmatmul
 ( int dim, double **Ahi, double **Alo, double **Fhi, double **Flo );
/*
 * DESCRIPTION :
 *   Replaces A with the product F*A.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in A and F;
 *   Ahi      high doubles of a matrix of dimension dim;
 *   Alo      low doubles of a matrix of dimension dim;
 *   Fhi      high doubles of a matrix of dimension dim;
 *   Flo      low doubles of a matrix of dimension dim.
 *
 * ON RETURN :
 *   A        high doubles of the product of F with A;
 *   A        low doubles of the product of F with A. */

void CPU_dbl2_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Uhi, double **Ulo,
   double *bhi, double *blo, double *xhi, double *xlo, double *lapsec );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Uhi      high doubles of an upper triangular matrix of dimension dim;
 *   Ulo      low doubles of an upper triangular matrix of dimension dim;
 *   bhi      high doubles of the right hand side of the linear system;
 *   blo      low doubles of the right hand side of the linear system;
 *   xhi      space allocated for dim doubles;
 *   xlo      space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhi      high doubles of the solution to the system U*x = b;
 *   xlo      low doubles of the solution to the system U*x = b;
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx2_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *lapsec );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Urehi    high doubles of the real parts of U;
 *   Urelo    low doubles of the real parts of U;
 *   Uimhi    high doubles of the imaginary parts of U;
 *   Uimlo    low doubles of the imaginary parts of U;
 *   brehi    high doubles of the real parts of b;
 *   brelo    low doubles of the real parts of b;
 *   bimhi    high doubles of the imaginary parts of b;
 *   bimlo    low doubles of the imaginary parts of b;
 *   xrehi    space allocated for dim doubles;
 *   xrelo    space allocated for dim doubles;
 *   ximhi    space allocated for dim doubles;
 *   ximlo    space allocated for dim doubles.
 *
 * ON RETURN :
 *   Urehi    the diagonal tiles contain the high doubles of
 *            the real parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Urelo    the diagonal tiles contain the low doubles of
 *            the real parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Uimhi    the diagonal tiles contain the high doubles of
 *            the imaginary parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   Uimlo    the diagonal tiles contain the low doubles of
 *            the imaginary parts of inverse matrices,
 *            for comparison with the result computed on the GPU;
 *   xrehi    high doubles of the real parts of the solution x;
 *   xrelo    low doubles of the real parts of the solution x;
 *   ximhi    high doubles of the imaginary parts of the solution x;
 *   ximlo    low doubles of the imaginary parts of the solution x;
 *   lapsec   elapsed time in seconds. */

#endif
