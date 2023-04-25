// The file dbl2_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in double double precision.

#ifndef __dbl2_bals_kernels_h__
#define __dbl2_bals_kernels_h__

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
   double *xhi, double *xlo,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   int vrblvl );
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
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qhi      high doubles of the Q in the QR of the Jacobian matrix;
 *   Qlo      low doubles of the Q in the QR of the Jacobian matrix;
 *   Rhi      high doubles of the R in the QR of the Jacobian matrix;
 *   Rlo      low doubles of the R in the QR of the Jacobian matrix;
 *   xhi      high doubles of the least squares solution;
 *   xlo      low doubles of the least squares solution;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions. */

void GPU_cmplx2_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo, 
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo,
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
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   vrblvl   is the verbose level, if zero, then no output.
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
 *   ximlo    low doubles of the imaginary parts of the solution;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions. */

void write_dbl2_qtbflops ( int ctype, int ncols, float lapsms );
/*
 * DESCRIPTION :
 *   Writes the number of flops and arithmetic intensity for Q^T*b.
 *
 * ON ENTRY :
 *   ctype    0 if real, 1 if complex;
 *   ncols    number of columns in the matrix;
 *   lapsms   time elapsed in milliseconds. */ 

void GPU_dbl2_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhi, double **Qlo, double *bhi, double *blo,
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
 *   Qhi      high doubles of the Q of the QR factorization;
 *   Qlo      low doubles of the Q of the QR factorization;
 *   bhi      high doubles of the right hand side vector of the system;
 *   blo      low doubles of the right hand side vector of the system;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   bhi      high doubles of the product of Q^T with b;
 *   blo      low doubles of the product of Q^T with b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs. */

void GPU_cmplx2_bals_qhb
 ( int ncols, int szt, int nbt,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *totqtblapsedms, int vrblvl );
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
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   brehi    high doubles of the real parts of Q^H*b;
 *   brelo    low doubles of the real parts of Q^H*b;
 *   bimhi    high doubles of the imaginary parts of Q^H*b;
 *   bimlo    low doubles of the imaginary parts of Q^H*b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs. */

void GPU_dbl2_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathi, double ***matlo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double **rhshi, double **rhslo,
   double **solhi, double **sollo,
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
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   vrblvl   the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qhi      high doubles of the Q of the QR of the Jacobian;
 *   Qlo      low doubles of the Q of the QR of the Jacobian;
 *   Rhi      high doubles of the R of the QR of the Jacobian;
 *   Rlo      low doubles of the R of the QR of the Jacobian;
 *   rhshi    high doubles of updated right hand side vectors;
 *   rhslo    low doubles of updated right hand side vectors;
 *   solhi    high doubles of the solution series;
 *   sollo    low doubles of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates. */

void GPU_cmplx2_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo, 
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
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   vrblvl   is the verbose level, if zero, then no output.
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
 *   solimlo  low doubles of the imaginary parts of the solution series;
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
