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

__global__ void cmplx_bals_tail
 ( int ncols, int szt, double *Are, double *Aim,
   double *xre, double *xim, double *bre, double *bim );
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
 *   Are      real parts of nrows-by-ncols matrix;
 *   Aim      imaginary parts of nrows-by-ncols matrix;
 *   xre      real parts of vector of dimension ncols;
 *   xim      imaginary parts of vector of dimension ncols;
 *   bre      real parts of vector of dimension nrows;
 *   bim      imaginary parts of vector of dimension nrows.
 *
 * ON RETURN :
 *   bre      real parts of updated right hand side vector;
 *   bim      imaginary parts of updated right hand side vector. */

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

__global__ void cmplx_bals_qhb
 ( int ncols, int szt, double *QHre, double *QHim,
   double *bre, double *bim, double *rre, double *rim );
/*
 * DESCRIPTION :
 *   Multiplies the Hermitian transpose of Q with b.
 *
 * REQUIRED : ncols = szt times the number of blocks,
 *   where ncols in the number of rows and columns in Qt 
 *   and the dimension of b.
 *
 * ON ENTRY :
 *   ncols    number of columns in Qt and the dimension of b;
 *   szt      size of each block (and tile);
 *   QHre     real parts of ncols-by-ncols matrix,
 *            the rows of QH contain the Hermitian transpose of Q;
 *   QHim     imaginary parts of ncols-by-ncols matrix,
 *            the rows of QH contain the Hermitian transpose of Q;
 *   bre      real parts of vector of dimension ncols;
 *   bim      imaginary parts of vector of dimension ncols;
 *   rre      real parts of vector of dimension ncols;
 *   rim      imaginary parts of vector of dimension ncols.
 *
 * ON RETURN :
 *   rre      real parts of product of QH with b;
 *   rim      imaginary parts of product of QH with b. */

void GPU_dbl_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, double *b, double *x, bool verbose );
/*
 * DESCRIPTION :
 *   Solves the head linear system in the least squares sense,
 *   with a QR factorization followed by a back substitution,
 *   wrapping the kernel launches for the blocked Householder QR
 *   followed by the tiled back substitution, on real data.
 *
 * REQUIRED : ncols = szt*nbt.
 *
 * ON ENTRY :
 *   nrows    number of rows in A and the dimension of b;
 *   ncols    number of columns in A and the dimension of x;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   A        an nrows-by-ncols matrix;
 *   Q        space allocated for a matrix of dimension nrows
 *   R        space allocated for a nrows-by-ncols matrix;
 *   b        defines the right hand side of the linear system;
 *   x        space for ncols numbers;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Q        the Q of the QR factorization of the Jacobian matrix;
 *   R        the R of the QR factorization of the Jacobian matrix;
 *   x        least squares solution. */

void GPU_cmplx_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Are, double **Aim, double **Qre, double **Qim,
   double **Rre, double **Rim, double *bre, double *bim,
   double *xre, double *xim, bool verbose );
/*
 * DESCRIPTION :
 *   Solves the head linear system in the least squares sense,
 *   with a QR factorization followed by a back substitution,
 *   wrapping the kernel launches for the blocked Householder QR
 *   followed by the tiled back substitution, on complex data.
 *
 * REQUIRED : ncols = szt*nbt.
 *
 * ON ENTRY :
 *   nrows    number of rows in A and the dimension of b;
 *   ncols    number of columns in A and the dimension of x;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   Are      real parts of an nrows-by-ncols matrix;
 *   Aim      real parts of an nrows-by-ncols matrix;
 *   Qre      space allocated for a matrix of dimension nrows
 *   Qim      space allocated for a matrix of dimension nrows
 *   Rre      space allocated for a nrows-by-ncols matrix;
 *   Rim      space allocated for a nrows-by-ncols matrix;
 *   bre      real parts of the right hand side of the linear system;
 *   bim      imaginary parts of the right hand side;
 *   xre      space for ncols numbers;
 *   xim      space for ncols numbers;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qre      real parts of the Q of the QR of the Jacobian;
 *   Qim      imaginary parts of the Q of the QR of the Jacobian;
 *   Rre      real parts of the R of the QR of the Jacobian;
 *   Rim      imaginary parts of the R of the QR of the Jacobian;
 *   xre      real parts of the least squares solution;
 *   xim      imaginary parts of the least squares solution. */

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

void GPU_cmplx_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matre, double ***matim, double **rhsre, double **rhsim,
   double **solre, double **solim, bool verbose );
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
 *   matre    real parts of the matrices of the linearized series;
 *   matim    imaginary parts of the matrices of the linearized series;
 *   rhsre    real part of the right hand sides;
 *   rhsim    imaginary parts of the right hand sides;
 *   solre    real parts of the solution computed up to stage-1;
 *   solim    imaginary parts of the solution computed up to stage-1;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   rhsre    real parts of the updated right hand sides;
 *   rhsim    imaginary parts of the updated right hand sides. */

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
 *   Q        the Q of the QR factorization;
 *   b        right hand side vector of the linear system;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   b        the product of Q^T with b. */

void GPU_cmplx_bals_qhb
 ( int ncols, int szt, int nbt, double **Qre, double **Qim,
   double *bre, double *bim, bool verbose );
/*
 * DESCRIPTION :
 *   The updated right hand side vector b is multiplied with Q^H.
 *
 * REQUIRED : ncols = szt*nbt.
 *
 * ON ENTRY :
 *   ncols    number of columns and rows in Q and the dimension
 *            of the vectors b and qtb;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   Qre      real parts of the Q of the QR factorization;
 *   Qim      imaginary parts the Q of the QR factorization;
 *   bre      real parts of the right hand side vector;
 *   bim      imaginary parts of the right hand side vector;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   b        real parts of the product of Q^H with b;
 *   b        imaginary parts of the product of Q^H with b. */

void GPU_dbl_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***mat, double **Q, double **R, double **rhs, double **sol,
   int vrblvl );
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
 *   Q        space for the Q of the QR factorization of the Jacobian;
 *   R        space for the R of the QR factorization of the Jacobian;
 *   rhs      degp1 vectors of dimension dim;
 *   sol      space allocated for degp1 vectors of dimension dim;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Q        the Q of the QR factorization of the Jacobian matrix;
 *   R        the R of the QR factorization of the Jacobian matrix;
 *   rhs      updated right hand side vectors;
 *   sol      coefficients of the solution series. */

void GPU_cmplx_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***matre, double ***matim, double **Qre, double **Qim,
   double **Rre, double **Rim, double **rhsre, double **rhsim,
   double **solre, double **solim, int vrblvl );
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
 *   matre    degp1 matrices of dimension dim;
 *   matim    degp1 matrices of dimension dim;
 *   Qre      space for the real parts of the Q of the QR;
 *   Qim      space for the imaginary parts of the Q of the QR;
 *   Rre      space for the real parts of the R of the QR;
 *   Rim      space for the imaginary parts of the R of the QR;
 *   rhsre    degp1 vectors of dimension dim;
 *   rhsim    degp1 vectors of dimension dim;
 *   solre    space allocated for degp1 vectors of dimension dim;
 *   solim    space allocated for degp1 vectors of dimension dim;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qre      real parts of the Q of the QR of the Jacobian matrix;
 *   Qim      imaginary parts of the Q of the QR of the Jacobian matrix;
 *   Rre      real parts of the R of the QR of the Jacobian matrix;
 *   Rim      imaginary parts of the R of the QR of the Jacobian matrix;
 *   rhsre    real parts of the updated right hand side vectors;
 *   rhsim    imaginary parts of the updated right hand side vectors;
 *   solre    real parts of the solution series;
 *   solim    imaginary parts of the solution series. */

#endif
