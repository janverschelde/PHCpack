// The file dbl2_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in double double precision.

#ifndef __dbl2_bals_kernels_h__
#define __dbl2_bals_kernels_h__

__global__ void dbl2_bals_tail
 ( int ncols, int szt, double *Ahi, double *Alo,
   double *xhi, double *xlo, double *bhi, double *blo );
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
 *   Ahi      nrows-by-ncols matrix of high doubles to multiply x with;
 *   Alo      nrows-by-ncols matrix of low doubles to multiply x with;
 *   xhi      high doubles of vector of dimension ncols;
 *   xlo      low doubles of vector of dimension ncols;
 *   bhi      high doubles of vector of dimension nrows;
 *   blo      low doubles of vector of dimension nrows.
 *
 * ON RETURN :
 *   bhi      high doubles of the updated right hand side vector;
 *   blo      low doubles of the updated right hand side vector. */

__global__ void dbl2_bals_qtb
 ( int ncols, int szt, double *Qthi, double *Qtlo,
   double *bhi, double *blo, double *rhi, double *rlo );
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
 *   Qthi     ncols-by-ncols matrix of high doubles to multiply x with,
 *            the rows of Qt contain the transpose of Q;
 *   Qtlo     ncols-by-ncols matrix of low doubles to multiply x with,
 *            the rows of Qt contain the transpose of Q;
 *   bhi      high doubles of vector of dimension ncols;
 *   blo      low doubles of vector of dimension ncols;
 *   rhi      high doubles of vector of dimension ncols.
 *   rlo      low doubles of vector of dimension ncols.
 *
 * ON RETURN :
 *   rhi      high doubles of product of Qt with b;
 *   rlo      low doubles of product of Qt with b. */

void GPU_dbl2_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahi, double **Alo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double *bhi, double *blo, 
   double *xhi, double *xlo, bool verbose );
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
 *   Qhi      space allocated for a matrix of dimension nrows
 *   Qlo      space allocated for a matrix of dimension nrows
 *   Rhi      space allocated for a nrows-by-ncols matrix;
 *   Rlo      space allocated for a nrows-by-ncols matrix;
 *   bhi      high doubles of the right hand side of the system;
 *   blo      low doubles of the right hand side of the system;
 *   xhi      space for ncols numbers;
 *   xlo      space for ncols numbers;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhi      high doubles of the Q in the QR of the Jacobian matrix;
 *   Qlo      low doubles of the Q in the QR of the Jacobian matrix;
 *   Rhi      high dobles of the  R in the QR f the Jacobian matrix;
 *   Rlo      low dobles of the  R in the QR f the Jacobian matrix;
 *   xhi      high dobles of the least squares solution;
 *   xlo      low doubles of the least squares solution. */

void GPU_dbl2_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathi, double ***matlo, double **rhshi, double **rhslo,
   double **solhi, double **sollo, bool verbose );
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
 *   mathi    high doubles of the matrices of the linearized series;
 *   matlo    low doubles of the matrices of the linearized series;
 *   rhshi    high doubles of the right hand sides of the system;
 *   rhslo    low doubles of the right hand sides of the system;
 *   solhi    high doubles of solution computed up to stage-1;
 *   sollo    low doubles of solution computed up to stage-1;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   rhshi    high doubles of updated right hand sides;
 *   rhslo    low doubles of updated right hand sides. */

void GPU_dbl2_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhi, double **Qlo, double *bhi, double *blo, bool verbose );
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
 *   Qhi      high doubles of the Q of the QR factorization;
 *   Qlo      low doubles of the Q of the QR factorization;
 *   bhi      high doubles of the right hand side vector of the system;
 *   blo      low doubles of the right hand side vector of the system;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   bhi      high doubles of the product of Q^T with b;
 *   blo      low doubles of the product of Q^T with b. */

void GPU_dbl2_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***mathi, double ***matlo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double **rhshi, double **rhslo,
   double **solhi, double **sollo, int vrblvl );
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
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   Qhi      space for the high doubles of the Q of the QR of the Jacobian;
 *   Qlo      space for the low doubles of the Q of the QR of the Jacobian;
 *   Rhi      space for the high doubles of the R of the QR fof the Jacobian;
 *   Rlo      space for the low doubles of the R of the QR fof the Jacobian;
 *   rhs      degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim;
 *   sollo    space allocated for degp1 vectors of dimension dim;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qhi      high doubles of the Q of the QR of the Jacobian;
 *   Qlo      low doubles of the Q of the QR of the Jacobian;
 *   Rhi      high doubles of the R of the QR of the Jacobian;
 *   Rlo      low doubles of the R of the QR of the Jacobian;
 *   rhshi    high doubles of updated right hand side vectors;
 *   rhslo    low doubles of updated right hand side vectors;
 *   solhi    high doubles of the solution series;
 *   sollo    low doubles of the solution series. */

#endif
