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

__global__ void cmplx2_bals_tail
 ( int ncols, int szt,
   double *Arehi, double *Arelo, double *Aimhi, double *Aimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo, 
   double *brehi, double *brelo, double *bimhi, double *bimlo );
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
 *   Arehi    high doubles of the real parts of nrows-by-ncols matrix;
 *   Arelo    low doubles of the real parts of nrows-by-ncols matrix;
 *   Aimhi    high doubles of the imaginary parts of nrows-by-ncols matrix;
 *   Aimlo    low doubles of the imaginary parts of nrows-by-ncols matrix;
 *   xrehi    high doubles of the real parts of vector of dimension ncols;
 *   xrelo    low doubles of the real parts of vector of dimension ncols;
 *   ximhi    high doubles of the imaginary parts of vector of dimension ncols;
 *   ximlo    low doubles of the imaginary parts of vector of dimension ncols;
 *   brehi    high doubles of the real parts of vector of dimension nrows;
 *   brelo    low doubles of the real parts of vector of dimension nrows;
 *   bimhi    high doubles of the imaginary parts of vector of dimension nrows;
 *   bimlo    low doubles of the imaginary parts of vector of dimension nrows.
 *
 * ON RETURN :
 *   brehi    high doubles of the real parts of updated right hand side;
 *   brelo    low doubles of the real parts of updated right hand side;
 *   bimhi    high doubles of the imag parts of updated right hand side;
 *   bimlo    low doubles of the imag parts of updated right hand side. */

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

__global__ void cmplx2_bals_qhb
 ( int ncols, int szt,
   double *QHrehi, double *QHrelo, double *QHimhi, double *QHimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *rrehi, double *rrelo, double *rimhi, double *rimlo );
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
 *   QHrehi   high doubles of the real parts of ncols-by-ncols matrix,
 *            the rows of QH contain the Hermitian transpose of Q;
 *   QHrelo   low doubles of the real parts of ncols-by-ncols matrix,
 *            the rows of QH contain the Hermitian transpose of Q;
 *   QHimhi   high doubles of the imaginary parts of ncols-by-ncols matrix,
 *            the rows of QH contain the Hermitian transpose of Q;
 *   QHimlo   low doubles of the imaginary parts of ncols-by-ncols matrix,
 *            the rows of QH contain the Hermitian transpose of Q;
 *   brehi    high doubles of the real parts of vector of dimension ncols;
 *   brelo    low doubles of the real parts of vector of dimension ncols;
 *   bimhi    high doubles of the imaginary parts of vector of dimension ncols;
 *   bimlo    low doubles of the imaginary parts of vector of dimension ncols;
 *   rrehi    high doubles of the real parts of vector of dimension ncols;
 *   rrelo    low doubles of the real parts of vector of dimension ncols;
 *   rimhi    high doubles of the imaginary parts of vector of dimension ncols.
 *   rimlo    low doubles of the imaginary parts of vector of dimension ncols.
 *
 * ON RETURN :
 *   rrehi    high doubles of the real parts of product of QH with b;
 *   rrelo    low doubles of the real parts of product of QH with b;
 *   rimhi    high doubles of the imaginary parts of product of QH with b;
 *   rimlo    low doubles of the imaginary parts of product of QH with b. */

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
 *   Rhi      high doubles of the R in the QR of the Jacobian matrix;
 *   Rlo      low doubles of the R in the QR of the Jacobian matrix;
 *   xhi      high doubles of the least squares solution;
 *   xlo      low doubles of the least squares solution. */

void GPU_cmplx2_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo, 
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   bool verbose );
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
 *   Arehi    high doubles of real parts of an nrows-by-ncols matrix;
 *   Arelo    low doubles of real parts of an nrows-by-ncols matrix;
 *   Aimhi    high doubles of real parts of an nrows-by-ncols matrix;
 *   Aimlo    low doubles of real parts of an nrows-by-ncols matrix;
 *   Qrehi    space allocated for a matrix of dimension nrows
 *   Qrelo    space allocated for a matrix of dimension nrows
 *   Qimhi    space allocated for a matrix of dimension nrows
 *   Qimlo    space allocated for a matrix of dimension nrows
 *   Rrehi    space allocated for a nrows-by-ncols matrix;
 *   Rrelo    space allocated for a nrows-by-ncols matrix;
 *   Rimhi    space allocated for a nrows-by-ncols matrix;
 *   Rimlo    space allocated for a nrows-by-ncols matrix;
 *   brehi    high doubles of real parts of the right hand side;
 *   brelo    low doubles of real parts of the right hand side;
 *   bimhi    high doubles of imaginary parts of the right hand side;
 *   bimlo    low doubles of imaginary parts of the right hand side;
 *   xrehi    space for ncols numbers;
 *   xrelo    space for ncols numbers;
 *   ximhi    space for ncols numbers;
 *   ximlo    space for ncols numbers;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qrehi    high doubles of the real parts of the Q of the QR;
 *   Qrelo    low doubles of the real parts of the Q of the QR;
 *   Qimhi    high doubles of the imaginary parts of the Q of the QR;
 *   Qimlo    low doubles of the imaginary parts of the Q of the QR;
 *   Rrehi    high doubles of the real parts of the R of the QR;
 *   Rrelo    low doubles of the real parts of the R of the QR;
 *   Rimhi    high doubles of the imaginary parts of the R of the QR;
 *   Rimlo    low doubles of the imaginary parts of the R of the QR;
 *   xrehi    high doubles of the real parts of the solution;
 *   xrelo    low doubles of the real parts of the solution;
 *   ximhi    high doubles of the imaginary parts of the solution;
 *   ximlo    low doubles of the imaginary parts of the solution. */

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

void GPU_cmplx2_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   bool verbose );
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
 *   matrehi  high doubles of the real parts of the matrices
 *            of the linearized series;
 *   matrelo  low doubles of the real parts of the matrices
 *            of the linearized series;
 *   matimhi  high doubles of the imaginary parts of the matrices
 *            of the linearized series;
 *   matimlo  low doubles of the imaginary parts of the matrices
 *            of the linearized series;
 *   rhsrehi  high doubles of the real part of the right hand sides;
 *   rhsrelo  low doubles of the real part of the right hand sides;
 *   rhsimhi  high doubles of the imaginary parts of the right hand sides;
 *   rhsimlo  low doubles of the imaginary parts of the right hand sides;
 *   solrehi  high doubles of the real parts of the solution,
 *            computed up to stage-1;
 *   solrelo  low doubles of the real parts of the solution,
 *            computed up to stage-1;
 *   solimhi  high doubles of the imaginary parts of the solution,
 *            computed up to stage-1;
 *   solimlo  low doubles of the imaginary parts of the solution,
 *            computed up to stage-1;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   rhsrehi  high doubles of the real parts of the updated right hand sides;
 *   rhsrelo  low doubles of the real parts of the updated right hand sides;
 *   rhsimhi  high doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimlo  low doubles of the imaginary parts
 *            of the updated right hand sides. */

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

void GPU_cmplx2_bals_qhb
 ( int ncols, int szt, int nbt,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo, bool verbose );
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
 *   Qrehi    high doubles of the real parts of the Q of the QR;
 *   Qrelo    low doubles of the real parts of the Q of the QR;
 *   Qimhi    high doubles of the imaginary parts the Q of the QR;
 *   Qimlo    low doubles of the imaginary parts the Q of the QR;
 *   brehi    high doubles of the real parts of the right hand side;
 *   brelo    low doubles of the real parts of the right hand side;
 *   bimhi    high doubles of the imaginary parts of the right hand side;
 *   bimlo    low doubles of the imaginary parts of the right hand side;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   brehi    high doubles of the real parts of Q^H*b;
 *   brelo    low doubles of the real parts of Q^H*b;
 *   bimhi    high doubles of the imaginary parts of Q^H*b;
 *   bimlo    low doubles of the imaginary parts of Q^H*b. */

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
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
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

void GPU_cmplx2_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo, 
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
 *   matrehi  degp1 matrices of dimension dim;
 *   matrelo  degp1 matrices of dimension dim;
 *   matimhi  degp1 matrices of dimension dim;
 *   matimlo  degp1 matrices of dimension dim;
 *   Qrehi    space for the high doubles of the real parts of the Q of the QR;
 *   Qrelo    space for the low doubles of the real parts of the Q of the QR;
 *   Qimhi    space for the high doubles of the imag parts of the Q of the QR;
 *   Qimlo    space for the low doubles of the imag parts of the Q of the QR;
 *   Rrehi    space for the high doubles of the real parts of the R of the QR;
 *   Rrelo    space for the low doubles of the real parts of the R of the QR;
 *   Rimhi    space for the high doubles of the imag parts of the R of the QR;
 *   Rimlo    space for the low doubles of the imag parts of the R of the QR;
 *   rhsrehi  degp1 vectors of dimension dim;
 *   rhsrelo  degp1 vectors of dimension dim;
 *   rhsimhi  degp1 vectors of dimension dim;
 *   rhsimlo  degp1 vectors of dimension dim;
 *   solrehi  space allocated for degp1 vectors of dimension dim;
 *   solrelo  space allocated for degp1 vectors of dimension dim;
 *   solimhi  space allocated for degp1 vectors of dimension dim;
 *   solimlo  space allocated for degp1 vectors of dimension dim;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qrehi    high doubles of the real parts of the Q of the QR;
 *   Qrelo    low doubles of the real parts of the Q of the QR;
 *   Qimhi    high doubles of the imaginary parts of the Q of the QR;
 *   Qimlo    low doubles of the imaginary parts of the Q of the QR;
 *   Rrehi    high doubles of the real parts of the R of the QR;
 *   Rrelo    low doubles of the real parts of the R of the QR;
 *   Rimhi    high doubles of the imaginary parts of the R of the QR;
 *   Rimlo    low doubles of the imaginary parts of the R of the QR;
 *   rhsrehi  high doubles of the real parts of the updated right hand side;
 *   rhsrelo  low doubles of the real parts of the updated right hand side;
 *   rhsimhi  high doubles of the imag parts of the updated right hand side;
 *   rhsimlo  low doubles of the imag parts of the updated right hand side;
 *   solrehi  high doubles of the real parts of the solution series;
 *   solrelo  low doubles of the real parts of the solution series;
 *   solimhi  high doubles of the imaginary parts of the solution series;
 *   solimlo  low doubles of the imaginary parts of the solution series. */

void GPU_dbl2_linear_residue
 ( int dim, int degp1, int szt, int nbt,
   double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear power series system.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 right hand side vectors of dimension dim;
 *   rhslo    degp1 right hand side vectors of dimension dim;
 *   solhi    degp1 solution vectors of dimension dim;
 *   sollo    degp1 solution vectors of dimension dim;
 *   resvechi has space for the residual power series;
 *   resveclo has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvechi are the high doubles of the residual power series;
 *   resveclo are the low doubles of the residual power series;
 *   resmaxhi is the high double of the maximum component in resvec;
 *   resmaxlo is the low double of the maximum component in resvec. */

void GPU_cmplx2_linear_residue
 ( int dim, int degp1, int szt, int nbt,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo, 
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear power series system.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   matrehi  degp1 matrices of dimension dim;
 *   matrelo  degp1 matrices of dimension dim;
 *   matimhi  degp1 matrices of dimension dim;
 *   matimlo  degp1 matrices of dimension dim;
 *   rhsrehi  degp1 right hand side vectors of dimension dim;
 *   rhsrelo  degp1 right hand side vectors of dimension dim;
 *   rhsimhi  degp1 right hand side vectors of dimension dim;
 *   rhsimlo  degp1 right hand side vectors of dimension dim;
 *   solrehi  degp1 solution vectors of dimension dim;
 *   solrelo  degp1 solution vectors of dimension dim;
 *   solimhi  degp1 solution vectors of dimension dim;
 *   solimlo  degp1 solution vectors of dimension dim;
 *   resvecrehi has space for the residual power series;
 *   resvecrelo has space for the residual power series;
 *   resvecimhi has space for the residual power series;
 *   resvecimlo has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvecrehi are the high doubles of the real parts of the residual series;
 *   resvecrelo are the low doubles of the real parts of the residual series;
 *   resvecimhi are the high doubles of the imag parts the residual series;
 *   resvecimlo are the low doubles of the imag parts the residual series;
 *   resmaxhi is the high double of the max norm of the residual series;
 *   resmaxlo is the low double of the max norm of the residual series. */

#endif
