// The file dbl_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in double precision.

#ifndef __dbl_bals_kernels_h__
#define __dbl_bals_kernels_h__

void GPU_dbl_qrbs_solve
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, double *b, double *x, bool verbose );
/*
 * DESCRIPTION :
 *   Solves the linear system in the least squares sense,
 *   with a QR factorization followed by a back substitution.
 *
 * REQUIRED : ncols = szt*nbt.
 *
 * ON ENTRY :
 *   nrows    number of rows in A and the dimension of b;
 *   ncols    number of columns in A and the dimension of x;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   Q        space allocated for a matrix of dimension nrows
 *   R        space allocated for a nrows-by-ncols matrix;
 *   b        defines the right hand side of the linear system;
 *   x        space for ncols numbers;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Q        the Q in a QR factorization of the Jacobian matrix;
 *   R        the R in a QR factorization of the Jacobian matrix;
 *   x        least squares solution. */

void GPU_dbl_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***mat, double **rhs, double **sol, int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format,
 *   using QR factorization and substitutions.
 *
 * REQUIRED : dim = szt*nbt.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   mat      degp1 matrices of dimension dim;
 *   rhs      degp1 vectors of dimension dim;
 *   sol      space allocated for degp1 vectors of dimension dim;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   sol      coefficients of the solution series. */

#endif
