// The file dbl4_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in quad double precision.

#ifndef __dbl4_bals_kernels_h__
#define __dbl4_bals_kernels_h__

__global__ void dbl4_bals_qtb
 ( int ncols, int szt,
   double *Qthihi, double *Qtlohi, double *Qthilo, double *Qtlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *rhihi, double *rlohi, double *rhilo, double *rlolo );
/*
 * DESCRIPTION :
 *   Multiplies the transpose of Q with b, on real data.
 *
 * REQUIRED : ncols = szt times the number of blocks,
 *   where ncols in the number of rows and columns in Qt 
 *   and the dimension of b.
 *
 * ON ENTRY :
 *   ncols    number of columns in Qt and the dimension of b;
 *   szt      size of each block (and tile);
 *   Qthihi   ncols-by-ncols matrix of highest doubles to multiply x with,
 *            the rows of Qt contain the transpose of Q;
 *   Qtlohi   ncols-by-ncols matrix of second highest doubles of Q;
 *   Qthilo   ncols-by-ncols matrix of second lowest doubles of Q;
 *   Qtlolo   ncols-by-ncols matrix of lowest doubles of Q;
 *   bhihi    highest doubles of vector of dimension ncols;
 *   blohi    second highest doubles of vector of dimension ncols;
 *   bhilo    second lowest doubles of vector of dimension ncols;
 *   blolo    lowest doubles of vector of dimension ncols;
 *   rhihi    highest doubles of vector of dimension ncols;
 *   rlohi    second highest doubles of vector of dimension ncols;
 *   rhilo    second lowest doubles of vector of dimension ncols;
 *   rlolo    lowest doubles of vector of dimension ncols.
 *
 * ON RETURN :
 *   rhihi    highest doubles of product of Qt with b;
 *   rlohi    second highest doubles of product of Qt with b;
 *   rhilo    second lowest doubles of product of Qt with b;
 *   rlolo    lowest doubles of product of Qt with b. */

__global__ void cmplx4_bals_qhb
 ( int ncols, int szt,
   double *QHrehihi, double *QHrelohi, double *QHrehilo, double *QHrelolo, 
   double *QHimhihi, double *QHimlohi, double *QHimhilo, double *QHimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *rrehihi, double *rrelohi, double *rrehilo, double *rrelolo,
   double *rimhihi, double *rimlohi, double *rimhilo, double *rimlolo );
/*
 * DESCRIPTION :
 *   Multiplies the Hermitian transpose of Q with b, on complex data.
 *
 * REQUIRED : ncols = szt times the number of blocks,
 *   where ncols in the number of rows and columns in Qt 
 *   and the dimension of b.
 *
 * ON ENTRY :
 *   ncols    number of columns in Qt and the dimension of b;
 *   szt      size of each block (and tile);
 *   QHrehihi are the highest doubles of the real parts of ncols-by-ncols
 *            matrix, the rows of QH contain the Hermitian transpose of Q;
 *   QHrelohi are the second highest doubles of the real parts of Q^H;
 *   QHrehilo are the second lowest doubles of the real parts of Q^H;
 *   QHrelolo are the lowest doubles of the real parts of Q^H;
 *   QHimhihi are the highest doubles of the imaginary parts of Q^H;
 *   QHimlohi are the second highest doubles of the imaginary parts of Q^H;
 *   QHimhilo are the second lowest doubles of the imaginary parts of Q^H;
 *   QHimlolo are the lowest doubles of the imaginary parts of Q;
 *   brehihi are the highest doubles of the real parts of b;
 *   brelohi are the second highest doubles of the real parts of b;
 *   brehilo are the second lowest doubles of the real parts of b;
 *   brelolo are the lowest doubles of the real parts of b;
 *   bimhihi are the highest doubles of the imaginary parts of b;
 *   bimlohi are the second highest doubles of the imaginary parts of b;
 *   bimhilo are the second lowest doubles of the imaginary parts of b;
 *   bimlolo are the lowest doubles of the imaginary parts of b;
 *   rrehihi has space for the highest doubles of the real parts;
 *   rrelohi has space for the second highest doubles of the real parts;
 *   rrehilo has space for the second lowest doubles of the real parts;
 *   rrelolo has space for the lowest doubles of the real parts;
 *   rimhihi has space for the highest doubles of the imaginary parts;
 *   rimlohi has space for the second highest doubles of the imaginary parts;
 *   rimhilo has space for the second lowest doubles of the imaginary parts;
 *   rimlolo has space for the lowest doubles of the imaginary parts.
 *
 * ON RETURN :
 *   rrehihi are the highest doubles of the real parts of QH*b;
 *   rrelohi are the second highest doubles of the real parts of QH*b;
 *   rrehilo are the second lowest doubles of the real parts of QH*b;
 *   rrelolo are the lowest doubles of the real parts of QH*b;
 *   rimhihi are the highest doubles of the imaginary parts of QH*b;
 *   rimlohi are the second highest doubles of the imaginary parts of QH*b;
 *   rimhilo are the second lowest doubles of the imaginary parts of QH*b;
 *   rimlolo are the lowest doubles of the imaginary parts of QH*b. */

void GPU_dbl4_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo, 
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
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
 *   Qhihi    space allocated for a matrix of dimension nrows
 *   Qlohi    space allocated for a matrix of dimension nrows
 *   Qhilo    space allocated for a matrix of dimension nrows
 *   Qlolo    space allocated for a matrix of dimension nrows
 *   Rhihi    space allocated for a nrows-by-ncols matrix;
 *   Rlohi    space allocated for a nrows-by-ncols matrix;
 *   Rhilo    space allocated for a nrows-by-ncols matrix;
 *   Rlolo    space allocated for a nrows-by-ncols matrix;
 *   bhihi    highest doubles of the right hand side of the system;
 *   blohi    second highest doubles of the right hand side of the system;
 *   bhilo    second lowest doubles of the right hand side of the system;
 *   blolo    lowest doubles of the right hand side of the system;
 *   xhihi    space for ncols numbers;
 *   xlohi    space for ncols numbers;
 *   xhilo    space for ncols numbers;
 *   xlolo    space for ncols numbers;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of the Q in the QR of the Jacobian matrix;
 *   Qlohi    second highest doubles of the Q in the QR of the Jacobian;
 *   Qhilo    second lowest doubles of the Q in the QR of the Jacobian;
 *   Qlolo    lowest doubles of the Q in the QR of the Jacobian;
 *   Rhihi    highest doubles of the R in the QR of the Jacobian;
 *   Rlohi    second highest doubles of the R in the QR of the Jacobian;
 *   Rhilo    second lowest doubles of the R in the QR of the Jacobian;
 *   Rlolo    lowest doubles of the R in the QR of the Jacobian;
 *   xhihi    highest doubles of the least squares solution;
 *   xlohi    second highest doubles of the least squares solution;
 *   xhilo    second lowest doubles of the least squares solution;
 *   xlolo    lowest doubles of the least squares solution;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions. */

void GPU_cmplx4_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo, 
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
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
 *   Arehihi  highest doubles of real parts of A;
 *   Arelohi  second highest doubles of real parts of A;
 *   Arehilo  second lowest doubles of real parts of A;
 *   Arelolo  lowest doubles of real parts of A;
 *   Aimhihi  highest doubles of real parts of A;
 *   Aimlohi  second highest doubles of real parts of A;
 *   Aimhilo  second lowest doubles of real parts of A;
 *   Aimlolo  lowest doubles of real parts of A;
 *   Qrehihi  space allocated for a matrix of dimension nrows
 *   Qrelohi  space allocated for a matrix of dimension nrows
 *   Qrehilo  space allocated for a matrix of dimension nrows
 *   Qrelolo  space allocated for a matrix of dimension nrows
 *   Qimhihi  space allocated for a matrix of dimension nrows
 *   Qimlohi  space allocated for a matrix of dimension nrows
 *   Qimhilo  space allocated for a matrix of dimension nrows
 *   Qimlolo  space allocated for a matrix of dimension nrows
 *   Rrehihi  space allocated for a nrows-by-ncols matrix;
 *   Rrelohi  space allocated for a nrows-by-ncols matrix;
 *   Rrehilo  space allocated for a nrows-by-ncols matrix;
 *   Rrelolo  space allocated for a nrows-by-ncols matrix;
 *   Rimhihi  space allocated for a nrows-by-ncols matrix;
 *   Rimlohi  space allocated for a nrows-by-ncols matrix;
 *   Rimhilo  space allocated for a nrows-by-ncols matrix;
 *   brehihi  highest doubles of real parts of the right hand side b;
 *   brelohi  second highest doubles of real parts of b;
 *   brehilo  second lowest doubles of real parts of b;
 *   brelolo  lowest doubles of real parts of b;
 *   bimhihi  highest doubles of imaginary parts of b;
 *   bimlohi  second highest doubles of imaginary parts of b;
 *   bimhilo  second lowest doubles of imaginary parts of b;
 *   bimlolo  lowest doubles of imaginary parts of b;
 *   xrehihi  space for ncols numbers;
 *   xrelohi  space for ncols numbers;
 *   xrehilo  space for ncols numbers;
 *   xrelolo  space for ncols numbers;
 *   ximhihi  space for ncols numbers;
 *   ximlohi  space for ncols numbers;
 *   ximhilo  space for ncols numbers;
 *   ximlolo  space for ncols numbers;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of the Q of the QR;
 *   Qrelohi  second highest doubles of the real parts of the Q of the QR;
 *   Qrehilo  second lowest doubles of the real parts of the Q of the QR;
 *   Qrelolo  lowest doubles of the real parts of the Q of the QR;
 *   Qimhihi  highest doubles of the imaginary parts of the Q of the QR;
 *   Qimlohi  second highest doubles of the imaginary parts of the Q of the QR;
 *   Qimhilo  second lowest doubles of the imaginary parts of the Q of the QR;
 *   Qimlolo  lowest doubles of the imaginary parts of the Q of the QR;
 *   Rrehihi  highest doubles of the real parts of the R of the QR;
 *   Rrelohi  second highest doubles of the real parts of the R of the QR;
 *   Rrehilo  second lowest doubles of the real parts of the R of the QR;
 *   Rrelolo  lowest doubles of the real parts of the R of the QR;
 *   Rimhihi  highest doubles of the imaginary parts of the R of the QR;
 *   Rimlohi  second highest doubles of the imaginary parts of the R of the QR;
 *   Rimhilo  second lowest doubles of the imaginary parts of the R of the QR;
 *   Rimlolo  lowest doubles of the imaginary parts of the R of the QR;
 *   xrehihi  highest doubles of the real parts of the solution;
 *   xrelohi  second highest doubles of the real parts of the solution;
 *   xrehilo  second lowest doubles of the real parts of the solution;
 *   xrelolo  lowest doubles of the real parts of the solution;
 *   ximhihi  highest doubles of the imaginary parts of the solution;
 *   ximlohi  second highest doubles of the imaginary parts of the solution;
 *   ximhilo  second lowest doubles of the imaginary parts of the solution;
 *   ximlolo  lowest doubles of the imaginary parts of the solution;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions. */

void write_dbl4_qtbflops ( int ctype, int ncols, float lapsms );
/*
 * DESCRIPTION :
 *   Writes the number of flops and arithmetic intensity for Q^T*b.
 *
 * ON ENTRY :
 *   ctype    0 if real, 1 if complex;
 *   ncols    number of columns in the matrix;
 *   lapsms   time elapsed in milliseconds. */ 

void GPU_dbl4_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
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
 *   Qhihi    highest doubles of the Q of the QR factorization;
 *   Qlohi    second highest doubles of the Q of the QR factorization;
 *   Qhilo    second lowest doubles of the Q of the QR factorization;
 *   Qlolo    lowest doubles of the Q of the QR factorization;
 *   bhihi    highest doubles of the right hand side vector;
 *   blohi    second highest doubles of the right hand side vector;
 *   bhilo    second lowest doubles of the right hand side vector;
 *   blolo    lowest doubles of the right hand side vector;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   bhihi    highest doubles of the product of Q^T with b;
 *   blohi    second highest doubles of the product of Q^T with b;
 *   bhilo    second lowest doubles of the product of Q^T with b;
 *   blolo    lowest doubles of the product of Q^T with b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs. */

void GPU_cmplx4_bals_qhb
 ( int ncols, int szt, int nbt,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
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
 *   Qrehihi  highest doubles of the real parts of the Q of the QR;
 *   Qrelohi  second highest doubles of the real parts of the Q of the QR;
 *   Qrehilo  second lowest doubles of the real parts of the Q of the QR;
 *   Qrelolo  lowest doubles of the real parts of the Q of the QR;
 *   Qimhihi  highest doubles of the imaginary parts the Q of the QR;
 *   Qimlohi  second highest doubles of the imaginary parts the Q of the QR;
 *   Qimhilo  second lowest doubles of the imaginary parts the Q of the QR;
 *   Qimlolo  lowest doubles of the imaginary parts the Q of the QR;
 *   brehihi  highest doubles of the real parts of the right hand side b;
 *   brelohi  second highest doubles of the real parts of b;
 *   brehilo  second lowest doubles of the real parts of b;
 *   brelolo  lowest doubles of the real parts of b;
 *   bimhihi  highest doubles of the imaginary parts of b;
 *   bimlohi  second highest doubles of the imaginary parts of b;
 *   bimhilo  second lowest doubles of the imaginary parts of b; 
 *   bimlolo  lowest doubles of the imaginary parts of b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   brehihi  highest doubles of the real parts of Q^H*b;
 *   brelohi  second highest doubles of the real parts of Q^H*b;
 *   brehilo  second lowest doubles of the real parts of Q^H*b;
 *   brelolo  lowest doubles of the real parts of Q^H*b;
 *   bimhihi  highest doubles of the imaginary parts of Q^H*b;
 *   bimlohi  second highest doubles of the imaginary parts of Q^H*b;
 *   bimhilo  second lowest doubles of the imaginary parts of Q^H*b;
 *   bimlolo  lowest doubles of the imaginary parts of Q^H*b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs. */

void GPU_dbl4_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
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
 *   mathihi  degp1 matrices of dimension dim;
 *   matlohi  degp1 matrices of dimension dim;
 *   mathilo  degp1 matrices of dimension dim;
 *   matlolo  degp1 matrices of dimension dim;
 *   Qhihi    space for the highest doubles of the Q of the QR;
 *   Qlohi    space for the second highest doubles of the Q of the QR;
 *   Qhilo    space for the second lowest doubles of the Q of the QR;
 *   Qlolo    space for the lowest doubles of the Q of the QR;
 *   Rhihi    space for the highest doubles of the R of the QR;
 *   Rlohi    space for the second highest doubles of the R of the QR;
 *   Rhilo    space for the second lowest doubles of the R of the QR;
 *   Rlolo    space for the lowest doubles of the R of the QR;
 *   rhshihi  degp1 vectors of dimension dim;
 *   rhslohi  degp1 vectors of dimension dim;
 *   rhshilo  degp1 vectors of dimension dim;
 *   rhslolo  degp1 vectors of dimension dim;
 *   solhihi  space allocated for degp1 vectors of dimension dim;
 *   sollohi  space allocated for degp1 vectors of dimension dim;
 *   solhilo  space allocated for degp1 vectors of dimension dim;
 *   sollolo  space allocated for degp1 vectors of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of the Q of the QR of the Jacobian;
 *   Qlohi    second highest doubles of the Q of the QR of the Jacobian;
 *   Qhilo    second lowest doubles of the Q of the QR of the Jacobian;
 *   Qlolo    lowest doubles of the Q of the QR of the Jacobian;
 *   Rhihi    highest doubles of the R of the QR of the Jacobian;
 *   Rlohi    second highest doubles of the R of the QR of the Jacobian;
 *   Rhilo    second lowest doubles of the R of the QR of the Jacobian;
 *   Rlolo    lowest doubles of the R of the QR of the Jacobian;
 *   rhshihi  highest doubles of updated right hand side vectors;
 *   rhslohi  second highest doubles of updated right hand side vectors;
 *   rhshilo  second lowest doubles of updated right hand side vectors;
 *   rhslolo  lowest doubles of updated right hand side vectors;
 *   solhihi  highest doubles of the solution series;
 *   sollohi  second highest doubles of the solution series.
 *   solhilo  second lowest doubles of the solution series;
 *   sollolo  lowest doubles of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent updates. */

void GPU_cmplx4_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo, 
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi, 
   double **solimhilo, double **solimlolo,
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
 *   matrehihi are degp1 matrices of dimension dim;
 *   matrelohi are degp1 matrices of dimension dim;
 *   matrehilo are degp1 matrices of dimension dim;
 *   matrelolo are degp1 matrices of dimension dim;
 *   matimhihi are degp1 matrices of dimension dim;
 *   matimlohi are degp1 matrices of dimension dim;
 *   matimhilo are degp1 matrices of dimension dim;
 *   matimlolo are degp1 matrices of dimension dim;
 *   Qrehihi  space for the highest doubles of
 *            the real parts of the Q of the QR;
 *   Qrelohi  space for the second highest doubles of
 *            the real parts of the Q of the QR;
 *   Qrehilo  space for the second lowest doubles of 
 *            the real parts of the Q of the QR;
 *   Qrelolo  space for the lowest doubles of 
 *            the real parts of the Q of the QR;
 *   Qimhihi  space for the highest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimlohi  space for the second highest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimhilo  space for the second lowest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimlolo  space for the lowest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Rrehihi  space for the highest doubles of
 *            the real parts of the R of the QR;
 *   Rrelohi  space for the second highest doubles of
 *            the real parts of the R of the QR;
 *   Rrehilo  space for the second lowest doubles of 
 *            the real parts of the R of the QR;
 *   Rrelolo  space for the lowest doubles of 
 *            the real parts of the R of the QR;
 *   Rimhihi  space for the highest doubles of
 *            the imaginary parts of the R of the QR;
 *   Rimlohi  space for the second highest doubles of
 *            the imaginary parts of the R of the QR;
 *   Rimhilo  space for the second lowest doubles of 
 *            the imaginary parts of the R of the QR;
 *   Rimlolo  space for the lowest doubles of 
 *            the imaginary parts of the R of the QR;
 *   rhsrehihi are degp1 vectors of dimension dim;
 *   rhsrelohi are degp1 vectors of dimension dim;
 *   rhsrehilo are degp1 vectors of dimension dim;
 *   rhsrelolo are degp1 vectors of dimension dim;
 *   rhsimhihi are degp1 vectors of dimension dim;
 *   rhsimlohi are degp1 vectors of dimension dim;
 *   rhsimhilo are degp1 vectors of dimension dim;
 *   rhsimlolo are degp1 vectors of dimension dim;
 *   solrehihi has space allocated for degp1 vectors of dimension dim;
 *   solrelohi has space allocated for degp1 vectors of dimension dim;
 *   solrehilo has space allocated for degp1 vectors of dimension dim;
 *   solrelolo has space allocated for degp1 vectors of dimension dim;
 *   solimhihi has space allocated for degp1 vectors of dimension dim;
 *   solimlohi has space allocated for degp1 vectors of dimension dim;
 *   solimhilo has space allocated for degp1 vectors of dimension dim;
 *   solimlolo has space allocated for degp1 vectors of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of the Q of the QR;
 *   Qrelohi  second highest doubles of the real parts of the Q of the QR;
 *   Qrehilo  second lowest doubles of the real parts of the Q of the QR;
 *   Qrelolo  lowest doubles of the real parts of the Q of the QR;
 *   Qimhihi  highest doubles of the imaginary parts of the Q of the QR;
 *   Qimlohi  second highest doubles of the imaginary parts of the Q of the QR;
 *   Qimhilo  second lowest doubles of the imaginary parts of the Q of the QR;
 *   Qimlolo  lowest doubles of the imaginary parts of the Q of the QR;
 *   Rrehihi  highest doubles of the real parts of the R of the QR;
 *   Rrelohi  second highest doubles of the real parts of the R of the QR;
 *   Rrehilo  second lowest doubles of the real parts of the R of the QR;
 *   Rrelolo  lowest doubles of the real parts of the R of the QR;
 *   Rimhihi  highest doubles of the imaginary parts of the R of the QR;
 *   Rimlohi  second highest doubles of the imaginary parts of the R of the QR;
 *   Rimhilo  second lowest doubles of the imaginary parts of the R of the QR;
 *   Rimlolo  lowest doubles of the imaginary parts of the R of the QR;
 *   rhsrehihi are the highest doubles of the real parts of the updated rhs;
 *   rhsrelohi are the 2nd highest doubles of the real parts of the updated rhs;
 *   rhsrehilo are the 2nd lowest doubles of the real parts of the updated rhs;
 *   rhsrelolo are the lowest doubles of the real parts of the updated rhs;
 *   rhsimhihi are the highest doubles of the imag parts of the updated rhs;
 *   rhsimlohi are the 2nd highest doubles of the imag parts of the updated rhs;
 *   rhsimhilo are the 2nd lowest doubles of the imag parts of the updated rhs;
 *   rhsimlolo are the lowest doubles of the imag parts of the updated rhs;
 *   solrehihi are the highest doubles of the real parts of the solution;
 *   solrelohi are the 2nd highest doubles of the real parts of the solution;
 *   solrehilo are the 2nd lowest doubles of the real parts of the solution;
 *   solrelolo are the lowest doubles of the real parts of the solution;
 *   solimhihi are the highest doubles of the imaginary parts of the solution;
 *   solimlohi are the 2nd highest doubles of the imag parts of the solution;
 *   solimhilo are the 2nd lowest doubles of the imag parts of the solution;
 *   solimlolo are the lowest doubles of the imag parts of the solution;
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
