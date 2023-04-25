// The file dbl_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in double precision.

#ifndef __dbl_bals_kernels_h__
#define __dbl_bals_kernels_h__

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
   double **A, double **Q, double **R, double *b, double *x,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   int vrblvl );
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
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Q        the Q of the QR factorization of the Jacobian matrix;
 *   R        the R of the QR factorization of the Jacobian matrix;
 *   x        least squares solution;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions. */

void GPU_cmplx_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Are, double **Aim, double **Qre, double **Qim,
   double **Rre, double **Rim, double *bre, double *bim,
   double *xre, double *xim,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   int vrblvl );
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
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qre      real parts of the Q of the QR of the Jacobian;
 *   Qim      imaginary parts of the Q of the QR of the Jacobian;
 *   Rre      real parts of the R of the QR of the Jacobian;
 *   Rim      imaginary parts of the R of the QR of the Jacobian;
 *   xre      real parts of the least squares solution;
 *   xim      imaginary parts of the least squares solution;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions. */

void write_dbl_qtbflops ( int ctype, int ncols, float lapsms );
/*
 * DESCRIPTION :
 *   Writes the number of flops and arithmetic intensity for Q^T*b.
 *
 * ON ENTRY :
 *   ctype    0 if real, 1 if complex;
 *   ncols    number of columns in the matrix;
 *   lapsms   time elapsed in milliseconds. */ 

void GPU_dbl_bals_qtb
 ( int ncols, int szt, int nbt, double **Q, double *b,
   double *totqtblapsedms, int vrblvl );
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
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   b        the product of Q^T with b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs. */

void GPU_cmplx_bals_qhb
 ( int ncols, int szt, int nbt, double **Qre, double **Qim,
   double *bre, double *bim, double *totqtblapsedms, int vrblvl );
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
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   bre      real parts of the product of Q^H with b;
 *   bim      imaginary parts of the product of Q^H with b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs. */

void GPU_dbl_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mat, double **Q, double **R, double **rhs, double **sol,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl );
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
 *   tailidx  the index of the start of the update in the tail;
 *   mat      degp1 matrices of dimension dim;
 *   Q        space for the Q of the QR factorization of the Jacobian;
 *   R        space for the R of the QR factorization of the Jacobian;
 *   rhs      degp1 vectors of dimension dim;
 *   sol      space allocated for degp1 vectors of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Q        the Q of the QR factorization of the Jacobian matrix;
 *   R        the R of the QR factorization of the Jacobian matrix;
 *   rhs      updated right hand side vectors;
 *   sol      coefficients of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates. */

void GPU_cmplx_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matre, double ***matim, double **Qre, double **Qim,
   double **Rre, double **Rim, double **rhsre, double **rhsim,
   double **solre, double **solim,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl );
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
 *   tailidx  the index of the start of the update in the tail;
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
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totbslapsedms accumulates the milliseconds spent on updates;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qre      real parts of the Q of the QR of the Jacobian matrix;
 *   Qim      imaginary parts of the Q of the QR of the Jacobian matrix;
 *   Rre      real parts of the R of the QR of the Jacobian matrix;
 *   Rim      imaginary parts of the R of the QR of the Jacobian matrix;
 *   rhsre    real parts of the updated right hand side vectors;
 *   rhsim    imaginary parts of the updated right hand side vectors;
 *   solre    real parts of the solution series;
 *   solim    imaginary parts of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates. */

#endif
