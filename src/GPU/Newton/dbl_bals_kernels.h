// The file dbl_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in double precision.

#ifndef __dbl_bals_kernels_h__
#define __dbl_bals_kernels_h__

__global__ void dbl_bals_tail
 ( int ncols, int szt, double *A, double *x, double *b );
/*
 * DESCRIPTION :
 *   Subtracts from the right hand side b the product of A with x.
 *
 * REQUIRED : nrows = szt times the number of blocks,
 *   where nrows in the number of rows in A and the dimension of b.
 *
 * ON ENTRY :
 *   ncols    number of columns in A and the dimension of x;
 *   szt      size of each block (and tile);
 *   A        nrows-by-ncols matrix to multiply x with;
 *   x        vector of dimension ncols;
 *   b        vector of dimension nrows.
 *
 * ON RETURN :
 *   b        updated right hand side vector. */

__global__ void dbl_bals_qtb
 ( int ncols, int szt, double *Qt, double *b, double *r );
/*
 * DESCRIPTION :
 *   Multiplies the transpose of Q with b.
 *
 * REQUIRED : ncols = szt times the number of blocks,
 *   where ncols in the number of rows and columns in Qt 
 *   and the dimension of b.
 *
 * ON ENTRY :
 *   ncols    number of columns in Qt and the dimension of b;
 *   szt      size of each block (and tile);
 *   Qt       ncols-by-ncols matrix to multiply x with,
 *            the rows of Qt contain the transpose of Q;
 *   b        vector of dimension ncols;
 *   r        vector of dimension ncols.
 *
 * ON RETURN :
 *   r        product of Qt with b. */

void GPU_dbl_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, double *b, double *x, bool verbose );
/*
 * DESCRIPTION :
 *   Solves the head linear system in the least squares sense,
 *   with a QR factorization followed by a back substitution,
 *   wrapping the kernel launches for the blocked Householder QR
 *   followed by the tiled back substitution.
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

void GPU_dbl_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mat, double **rhs, double **sol, bool verbose );
/*
 * DESCRIPTION :
 *   After each block of coefficients of the series,
 *   kernels are launched for the multiplication of the tail matrices
 *   with the solution coefficients to update the right hand sides
 *   of the linear system of power series.
 *
 * REQUIRED : ncols = szt*nbt.
 *
 * ON ENTRY :
 *   nrows    number of rows in A and the dimension of b;
 *   ncols    number of columns in A and the dimension of x;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   degp1    degree plus one, total number of coefficient blocks;
 *   stage    coefficient blocks up to stage-1 are computed;
 *   mat      matrices of the linearized power series;
 *   rhs      right hand side vectors of the linear system;
 *   sol      solution coefficients computed up to stage-1;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   rhs      updated right hand sides. */

void GPU_dbl_bals_qtb
 ( int ncols, int szt, int nbt, double **Q, double *b, bool verbose );
/*
 * DESCRIPTION :
 *   The updated right hand side vector b is multiplied with Q^T.
 *
 * REQUIRED : ncols = szt*nbt.
 *
 * ON ENTRY :
 *   ncols    number of columns and rows in Q and the dimension
 *            of the vectors b and qtb;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   b        right hand side vector of the linear system;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   b        the product of Q^T with b. */

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
