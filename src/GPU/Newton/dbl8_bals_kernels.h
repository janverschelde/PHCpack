// The file dbl8_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in octo double precision.

#ifndef __dbl8_bals_kernels_h__
#define __dbl8_bals_kernels_h__

__global__ void dbl8_bals_qtb
 ( int ncols, int szt,
   double *Qthihihi, double *Qtlohihi, double *Qthilohi, double *Qtlolohi,
   double *Qthihilo, double *Qtlohilo, double *Qthilolo, double *Qtlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *rhihihi, double *rlohihi, double *rhilohi, double *rlolohi,
   double *rhihilo, double *rlohilo, double *rhilolo, double *rlololo );
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
 *   Qthihihi is ncols-by-ncols matrix of highest doubles of Q,
 *            the rows of Qt contain the transpose of Q;
 *   Qtlohihi is ncols-by-ncols matrix of second highest doubles of Q;
 *   Qthilohi is ncols-by-ncols matrix of third highest doubles of Q;
 *   Qtlolohi is ncols-by-ncols matrix of fourth highest doubles of Q;
 *   Qthihilo is ncols-by-ncols matrix of fourth lowest doubles of Q;
 *   Qtlohilo is ncols-by-ncols matrix of second highest doubles of Q;
 *   Qthilolo is ncols-by-ncols matrix of second lowest doubles of Q;
 *   Qtlololo is ncols-by-ncols matrix of lowest doubles of Q;
 *   bhihihi  highest doubles of vector of dimension ncols;
 *   blohihi  second highest doubles of vector of dimension ncols;
 *   bhilohi  third highest doubles of vector of dimension ncols;
 *   blolohi  fourth highest doubles of vector of dimension ncols;
 *   bhihilo  fourth lowest doubles of vector of dimension ncols;
 *   blohilo  third lowest doubles of vector of dimension ncols;
 *   bhilolo  second lowest doubles of vector of dimension ncols;
 *   blololo  lowest doubles of vector of dimension ncols;
 *   rhihihi  highest doubles of vector of dimension ncols;
 *   rlohihi  second highest doubles of vector of dimension ncols;
 *   rhilohi  third highest doubles of vector of dimension ncols;
 *   rlolohi  fourth highest doubles of vector of dimension ncols.
 *   rhihilo  fourth lowest doubles of vector of dimension ncols;
 *   rlohilo  third lowest doubles of vector of dimension ncols;
 *   rhilolo  second lowest doubles of vector of dimension ncols;
 *   rlololo  lowest doubles of vector of dimension ncols.
 *
 * ON RETURN :
 *   rhihihi  highest doubles of product of Qt with b;
 *   rlohihi  second highest doubles of product of Qt with b;
 *   rhilohi  third highest doubles of product of Qt with b;
 *   rlolohi  fourth highest doubles of product of Qt with b;
 *   rhihilo  fourth lowest doubles of product of Qt with b;
 *   rlohilo  third lowest doubles of product of Qt with b;
 *   rhilolo  second lowest doubles of product of Qt with b;
 *   rlololo  lowest doubles of product of Qt with b. */

__global__ void cmplx8_bals_qhb
 ( int ncols, int szt,
   double *QHrehihihi, double *QHrelohihi,
   double *QHrehilohi, double *QHrelolohi, 
   double *QHrehihilo, double *QHrelohilo,
   double *QHrehilolo, double *QHrelololo, 
   double *QHimhihihi, double *QHimlohihi,
   double *QHimhilohi, double *QHimlolohi,
   double *QHimhihilo, double *QHimlohilo,
   double *QHimhilolo, double *QHimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *rrehihihi, double *rrelohihi, double *rrehilohi, double *rrelolohi,
   double *rrehihilo, double *rrelohilo, double *rrehilolo, double *rrelololo,
   double *rimhihihi, double *rimlohihi, double *rimhilohi, double *rimlolohi,
   double *rimhihilo, double *rimlohilo, double *rimhilolo, double *rimlololo );
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
 *   QHrehihihi are the highest doubles of the real parts of ncols-by-ncols
 *            matrix, the rows of QH contain the Hermitian transpose of Q;
 *   QHrelohihi are the second highest doubles of the real parts of Q^H;
 *   QHrehilohi are the third highest doubles of the real parts of Q^H;
 *   QHrelolohi are the fourth highest doubles of the real parts of Q^H;
 *   QHrehihilo are the fourth lowest doubles of the real parts of Q^H;
 *   QHrelohilo are the third lowest doubles of the real parts of Q^H;
 *   QHrehilolo are the second lowest doubles of the real parts of Q^H;
 *   QHrelololo are the lowest doubles of the real parts of Q^H;
 *   QHimhihihi are the highest doubles of the imaginary parts of Q^H;
 *   QHimlohihi are the second highest doubles of the imaginary parts of Q^H;
 *   QHimhilohi are the third highest doubles of the imaginary parts of Q^H;
 *   QHimlolohi are the fourth highest doubles of the imaginary parts of Q^H;
 *   QHimhihilo are the fourth lowest doubles of the imaginary parts of Q^H;
 *   QHimlohilo are the third lowest doubles of the imaginary parts of Q^H;
 *   QHimhilolo are the second lowest doubles of the imaginary parts of Q^H;
 *   QHimlololo are the lowest doubles of the imaginary parts of Q;
 *   brehihihi are the highest doubles of the real parts of b;
 *   brelohihi are the second highest doubles of the real parts of b;
 *   brehilohi are the third highest doubles of the real parts of b;
 *   brelolohi are the fourth highest doubles of the real parts of b;
 *   brehihilo are the fourth lowest doubles of the real parts of b;
 *   brelohilo are the third lowest doubles of the real parts of b;
 *   brehilolo are the second lowest doubles of the real parts of b;
 *   brelololo are the lowest doubles of the real parts of b;
 *   bimhihihi are the highest doubles of the imaginary parts of b;
 *   bimlohihi are the second highest doubles of the imaginary parts of b;
 *   bimhilohi are the third highest doubles of the imaginary parts of b;
 *   bimlolohi are the fourth highest doubles of the imaginary parts of b;
 *   bimhihilo are the fourth lowest doubles of the imaginary parts of b;
 *   bimlohilo are the third lowest doubles of the imaginary parts of b;
 *   bimhilolo are the second lowest doubles of the imaginary parts of b;
 *   bimlololo are the lowest doubles of the imaginary parts of b;
 *   rrehihihi has space for the highest doubles of the real parts;
 *   rrelohihi has space for the second highest doubles of the real parts;
 *   rrehilohi has space for the third highest doubles of the real parts;
 *   rrelolohi has space for the fourth highest doubles of the real parts;
 *   rrehihilo has space for the fourth lowest doubles of the real parts;
 *   rrelohilo has space for the third lowest doubles of the real parts;
 *   rrehilolo has space for the second lowest doubles of the real parts;
 *   rrelololo has space for the lowest doubles of the real parts;
 *   rimhihihi has space for the highest doubles of the imaginary parts;
 *   rimlohihi has space for the second highest doubles of the imag parts;
 *   rimhilohi has space for the third highest doubles of the imag parts;
 *   rimlolohi has space for the fourth highest doubles of the imag parts;
 *   rimhihilo has space for the fourth lowest doubles of the imag parts;
 *   rimlohilo has space for the third lowest doubles of the imag parts;
 *   rimhilolo has space for the second lowest doubles of the imag parts;
 *   rimlololo has space for the lowest doubles of the imaginary parts.
 *
 * ON RETURN :
 *   rrehihihi are the highest doubles of the real parts of QH*b;
 *   rrelohihi are the second highest doubles of the real parts of QH*b;
 *   rrehilohi are the third highest doubles of the real parts of QH*b;
 *   rrelolohi are the fourth highest doubles of the real parts of QH*b;
 *   rrehihilo are the fourth lowest doubles of the real parts of QH*b;
 *   rrelohilo are the third lowest doubles of the real parts of QH*b;
 *   rrehilolo are the second lowest doubles of the real parts of QH*b;
 *   rrelololo are the lowest doubles of the real parts of QH*b;
 *   rimhihihi are the highest doubles of the imaginary parts of QH*b;
 *   rimlohihi are the second highest doubles of the imaginary parts of QH*b;
 *   rimhilohi are the third highest doubles of the imaginary parts of QH*b;
 *   rimlolohi are the fourth highest doubles of the imaginary parts of QH*b;
 *   rimhihilo are the fourth lowest doubles of the imaginary parts of QH*b;
 *   rimlohilo are the third lowest doubles of the imaginary parts of QH*b;
 *   rimhilolo are the second lowest doubles of the imaginary parts of QH*b;
 *   rimlololo are the lowest doubles of the imaginary parts of QH*b. */

void GPU_dbl8_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
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
 *   Qhihihi  space allocated for a matrix of dimension nrows
 *   Qlohihi  space allocated for a matrix of dimension nrows
 *   Qhilohi  space allocated for a matrix of dimension nrows
 *   Qlolohi  space allocated for a matrix of dimension nrows
 *   Qhihilo  space allocated for a matrix of dimension nrows
 *   Qlohilo  space allocated for a matrix of dimension nrows
 *   Qhilolo  space allocated for a matrix of dimension nrows
 *   Qlololo  space allocated for a matrix of dimension nrows
 *   Rhihihi  space allocated for a nrows-by-ncols matrix;
 *   Rlohihi  space allocated for a nrows-by-ncols matrix;
 *   Rhilohi  space allocated for a nrows-by-ncols matrix;
 *   Rlolohi  space allocated for a nrows-by-ncols matrix;
 *   Rhihilo  space allocated for a nrows-by-ncols matrix;
 *   Rlohilo  space allocated for a nrows-by-ncols matrix;
 *   Rhilolo  pace allocated for a nrows-by-ncols matrix;
 *   Rlololo  space allocated for a nrows-by-ncols matrix;
 *   bhihihi  highest doubles of the right hand side;
 *   blohihi  second highest doubles of the right hand side;
 *   bhilohi  third highest doubles of the right hand side;
 *   blolohi  fourth highest doubles of the right hand side;
 *   bhihilo  fourth lowest doubles of the right hand side;
 *   blohilo  third lowest doubles of the right hand side;
 *   bhilolo  second lowest doubles of the right hand side;
 *   blololo  lowest doubles of the right hand side;
 *   xhihihi  space for ncols numbers;
 *   xlohihi  space for ncols numbers;
 *   xhilohi  space for ncols numbers;
 *   xlolohi  space for ncols numbers;
 *   xhihilo  space for ncols numbers;
 *   xlohilo  space for ncols numbers;
 *   xhilolo  space for ncols numbers;
 *   xlololo  space for ncols numbers;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qhihihi  highest doubles of the Q in the QR of the Jacobian;
 *   Qlohihi  second highest doubles of the Q in the QR of the Jacobian;
 *   Qhilohi  third highest doubles of the Q in the QR of the Jacobian;
 *   Qlolohi  fourth highest doubles of the Q in the QR of the Jacobian;
 *   Qhihilo  fourth lowest doubles of the Q in the QR of the Jacobian;
 *   Qlohilo  third lowest doubles of the Q in the QR of the Jacobian;
 *   Qhilolo  second lowest doubles of the Q in the QR of the Jacobian;
 *   Qlololo  lowest doubles of the Q in the QR of the Jacobian;
 *   Rhihihi  highest doubles of the R in the QR of the Jacobian;
 *   Rlohihi  second highest doubles of the R in the QR of the Jacobian;
 *   Rhilohi  third highest doubles of the R in the QR of the Jacobian;
 *   Rlolohi  fourth highest doubles of the R in the QR of the Jacobian;
 *   Rhihilo  fourth lowest doubles of the R in the QR of the Jacobian;
 *   Rlohilo  third lowest doubles of the R in the QR of the Jacobian;
 *   Rhilolo  second lowest doubles of the R in the QR of the Jacobian;
 *   Rlololo  lowest doubles of the R in the QR of the Jacobian;
 *   xhihihi  highest doubles of the least squares solution;
 *   xlohihi  second highest doubles of the least squares solution;
 *   xhilohi  third highest doubles of the least squares solution;
 *   xlolohi  fourth highest doubles of the least squares solution;
 *   xhihilo  fourth lowest doubles of the least squares solution;
 *   xlohilo  third lowest doubles of the least squares solution;
 *   xhilolo  second lowest doubles of the least squares solution;
 *   xlololo  lowest doubles of the least squares solution;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions. */

void GPU_cmplx8_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi, 
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo, 
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
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
 *   Arehihihi are the highest doubles of real parts of A;
 *   Arelohihi are the second highest doubles of real parts of A;
 *   Arehilohi are the third highest doubles of real parts of A;
 *   Arelolohi are the fourth highest doubles of real parts of A;
 *   Arehihilo are the fourth lowest doubles of real parts of A;
 *   Arelohilo are the third lowest doubles of real parts of A;
 *   Arehilolo are the second lowest doubles of real parts of A;
 *   Arelololo are the lowest doubles of real parts of A;
 *   Aimhihihi are the highest doubles of real parts of A;
 *   Aimlohihi are the second highest doubles of real parts of A;
 *   Aimhilohi are the third highest doubles of real parts of A;
 *   Aimlolohi are the fourth highest doubles of real parts of A;
 *   Aimhihilo are the fourth lowest doubles of real parts of A;
 *   Aimlohilo are the third lowest doubles of real parts of A;
 *   Aimhilolo are the second lowest doubles of real parts of A;
 *   Aimlololo are the lowest doubles of real parts of A;
 *   Qrehihihi has space allocated for a matrix of dimension nrows
 *   Qrelohihi has space allocated for a matrix of dimension nrows
 *   Qrehilohi has space allocated for a matrix of dimension nrows
 *   Qrelolohi has space allocated for a matrix of dimension nrows
 *   Qrehihilo has space allocated for a matrix of dimension nrows
 *   Qrelohilo has space allocated for a matrix of dimension nrows
 *   Qrehilolo has space allocated for a matrix of dimension nrows
 *   Qrelololo has space allocated for a matrix of dimension nrows
 *   Qimhihihi has space allocated for a matrix of dimension nrows
 *   Qimlohihi has space allocated for a matrix of dimension nrows
 *   Qimhilohi has space allocated for a matrix of dimension nrows
 *   Qimlolohi has space allocated for a matrix of dimension nrows
 *   Qimhihilo has space allocated for a matrix of dimension nrows
 *   Qimlohilo has space allocated for a matrix of dimension nrows
 *   Qimhilolo has space allocated for a matrix of dimension nrows
 *   Qimlololo has space allocated for a matrix of dimension nrows
 *   Rrehihihi has space allocated for a nrows-by-ncols matrix;
 *   Rrelohihi has space allocated for a nrows-by-ncols matrix;
 *   Rrehilohi has space allocated for a nrows-by-ncols matrix;
 *   Rrelolohi has space allocated for a nrows-by-ncols matrix;
 *   Rrehihilo has space allocated for a nrows-by-ncols matrix;
 *   Rrelohilo has space allocated for a nrows-by-ncols matrix;
 *   Rrehilolo has space allocated for a nrows-by-ncols matrix;
 *   Rrelololo has space allocated for a nrows-by-ncols matrix;
 *   Rimhihihi has space allocated for a nrows-by-ncols matrix;
 *   Rimlohihi has space allocated for a nrows-by-ncols matrix;
 *   Rimhilohi has space allocated for a nrows-by-ncols matrix;
 *   Rimlolohi has space allocated for a nrows-by-ncols matrix;
 *   Rimhihilo has space allocated for a nrows-by-ncols matrix;
 *   Rimlohilo has space allocated for a nrows-by-ncols matrix;
 *   Rimhilolo has space allocated for a nrows-by-ncols matrix;
 *   Rimlololo has space allocated for a nrows-by-ncols matrix;
 *   brehihihi are the highest doubles of real parts of the right hand side b;
 *   brelohihi are the second highest doubles of real parts of b;
 *   brehilohi are the third highest doubles of real parts of b;
 *   brelolohi are the fourth highest doubles of real parts of b;
 *   brehihilo are the fourth lowest doubles of real parts of b;
 *   brelohilo are the third lowest doubles of real parts of b;
 *   brehilolo are the second lowest doubles of real parts of b;
 *   brelololo are the lowest doubles of real parts of b;
 *   bimhihihi are the highest doubles of imaginary parts of b;
 *   bimlohihi are the second highest doubles of imaginary parts of b;
 *   bimhilohi are the third highest doubles of imaginary parts of b;
 *   bimlolohi are the fourth highest doubles of imaginary parts of b;
 *   bimhihilo are the fourth lowest doubles of imaginary parts of b;
 *   bimlohilo are the third lowest doubles of imaginary parts of b;
 *   bimhilolo are the second lowest doubles of imaginary parts of b;
 *   bimlololo are the lowest doubles of imaginary parts of b;
 *   xrehihihi has space for ncols numbers;
 *   xrelohihi has space for ncols numbers;
 *   xrehilohi has space for ncols numbers;
 *   xrelolohi has space for ncols numbers;
 *   xrehihilo has space for ncols numbers;
 *   xrelohilo has space for ncols numbers;
 *   xrehilolo has space for ncols numbers;
 *   xrelololo has space for ncols numbers;
 *   ximhihihi has space for ncols numbers;
 *   ximlohihi has space for ncols numbers;
 *   ximhilohi has space for ncols numbers;
 *   ximlolohi has space for ncols numbers;
 *   ximhihilo has space for ncols numbers;
 *   ximlohilo has space for ncols numbers;
 *   ximhilolo has space for ncols numbers;
 *   ximlololo has space for ncols numbers;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qrehihihi are the highest doubles of the real parts of the Q of the QR;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third highest doubles of the real parts of Q;
 *   Qrelolohi are the fourth highest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo are the third lowest doubles of the real parts of Q;
 *   Qrehilolo are the second lowest doubles of the real parts of Q;
 *   Qrelololo are the lowest doubles of the real parts of Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi are the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q;
 *   Rrehihihi are the highest doubles of the real parts of the R of the QR;
 *   Rrelohihi are the second highest doubles of the real parts of R;
 *   Rrehilohi are the third highest doubles of the real parts of R;
 *   Rrelolohi are the fourth highest doubles of the real parts of R;
 *   Rrehihilo are the fourth lowest doubles of the real parts of R;
 *   Rrelohilo are the third lowest doubles of the real parts of R;
 *   Rrehilolo are the second lowest doubles of the real parts of R;
 *   Rrelololo are the lowest doubles of the real parts of R;
 *   Rimhihihi are the highest doubles of the imaginary parts of R;
 *   Rimlohihi are the second highest doubles of the imaginary parts of R;
 *   Rimhilohi are the third highest doubles of the imaginary parts of R;
 *   Rimlolohi are the fourth highest doubles of the imaginary parts of R;
 *   Rimhihilo are the fourth lowest doubles of the imaginary parts of R;
 *   Rimlohilo are the third lowest doubles of the imaginary parts of R;
 *   Rimhilolo are the second lowest doubles of the imaginary parts of R;
 *   Rimlololo are the lowest doubles of the imaginary parts of R;
 *   xrehihihi are the highest doubles of the real parts of the solution x;
 *   xrelohihi are the second highest doubles of the real parts of x;
 *   xrehilohi are the third highest doubles of the real parts of x;
 *   xrelolohi are the fourth highest doubles of the real parts of x;
 *   xrehihilo are the fourth lowest doubles of the real parts of x;
 *   xrelohilo are the third lowest doubles of the real parts of x;
 *   xrehilolo are the second lowest doubles of the real parts of x;
 *   xrelololo are the lowest doubles of the real parts of x;
 *   ximhihihi are the highest doubles of the imaginary parts of x;
 *   ximlohihi are the second highest doubles of the imaginary parts of x;
 *   ximhilohi are the third highest doubles of the imaginary parts of x;
 *   ximlolohi are the fourth highest doubles of the imaginary parts of x;
 *   ximhihilo are the fourth lowest doubles of the imaginary parts of x;
 *   ximlohilo are the third lowest doubles of the imaginary parts of x;
 *   ximhilolo are the second lowest doubles of the imaginary parts of x;
 *   ximlololo are the lowest doubles of the imaginary parts of x;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions. */

void write_dbl8_qtbflops ( int ctype, int ncols, float lapsms );
/*
 * DESCRIPTION :
 *   Writes the number of flops and arithmetic intensity for Q^T*b.
 *
 * ON ENTRY :
 *   ctype    0 if real, 1 if complex;
 *   ncols    number of columns in the matrix;
 *   lapsms   time elapsed in milliseconds. */ 

void GPU_dbl8_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *totqtblapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   The updated right hand side vector b is multiplied with Q^T,
 *   on real data.
 *
 * REQUIRED : ncols = szt*nbt.
 *
 * ON ENTRY :
 *   ncols    number of columns and rows in Q and the dimension
 *            of the vectors b and qtb;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   Qhihihi  highest doubles of the Q of the QR factorization;
 *   Qlohihi  second highest doubles of the Q of the QR factorization;
 *   Qhilohi  third highest doubles of the Q of the QR factorization;
 *   Qlolohi  fourth highest doubles of the Q of the QR factorization;
 *   Qhihilo  fourth lowest doubles of the Q of the QR factorization;
 *   Qlohilo  third lowest doubles of the Q of the QR factorization;
 *   Qhilolo  second lowest doubles of the Q of the QR factorization;
 *   Qlololo  lowest doubles of the Q of the QR factorization;
 *   bhihihi  highest doubles of the right hand side vector;
 *   blohihi  second highest doubles of the right hand side vector;
 *   bhilohi  third highest doubles of the right hand side vector;
 *   blolohi  fourth highest doubles of the right hand side vector;
 *   bhihilo  fourth lowest doubles of the right hand side vector;
 *   blohilo  third lowest doubles of the right hand side vector;
 *   bhilolo  second lowest doubles of the right hand side vector;
 *   blololo  lowest doubles of the right hand side vector;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   bhihihi  highest doubles of the product of Q^T with b;
 *   blohihi  second highest doubles of the product of Q^T with b;
 *   bhilohi  third highest doubles of the product of Q^T with b;
 *   blolohi  fourth highest doubles of the product of Q^T with b;
 *   bhihilo  fourth lowest doubles of the product of Q^T with b;
 *   blohilo  third lowest doubles of the product of Q^T with b;
 *   bhilolo  second lowest doubles of the product of Q^T with b;
 *   blololo  lowest doubles of the product of Q^T with b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs. */

void GPU_cmplx8_bals_qhb
 ( int ncols, int szt, int nbt,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *totqtblapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   The updated right hand side vector b is multiplied with Q^H,
 *   on complex data.
 *
 * REQUIRED : ncols = szt*nbt.
 *
 * ON ENTRY :
 *   ncols    number of columns and rows in Q and the dimension
 *            of the vectors b and qtb;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   Qrehihihi are the highest doubles of the real parts of the Q of the QR;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third highest doubles of the real parts of Q;
 *   Qrelolohi are the fourth highest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo are the third lowest doubles of the real parts of Q;
 *   Qrehilolo are the second lowest doubles of the real parts of Q;
 *   Qrelololo are the lowest doubles of the real parts of Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi are the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q;
 *   brehihihi are the highest doubles of the real parts of b;
 *   brelohihi are the second highest doubles of the real parts of b;
 *   brehilohi are the third highest doubles of the real parts of b;
 *   brelolohi are the fourth highest doubles of the real parts of b;
 *   brehihilo are the fourth lowest doubles of the real parts of b;
 *   brelohilo are the third lowest doubles of the real parts of b;
 *   brehilolo are the second lowest doubles of the real parts of b;
 *   brelololo are the lowest doubles of the real parts of b;
 *   bimhihihi are the highest doubles of the imaginary parts of b;
 *   bimlohihi are the second highest doubles of the imaginary parts of b;
 *   bimhilohi are the third highest doubles of the imaginary parts of b;
 *   bimlolohi are the fourth highest doubles of the imaginary parts of b;
 *   bimhihilo are the fourth lowest doubles of the imaginary parts of b; 
 *   bimlohilo are the third lowest doubles of the imaginary parts of b; 
 *   bimhilolo are the second lowest doubles of the imaginary parts of b; 
 *   bimlololo are the lowest doubles of the imaginary parts of b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   brehihihi are the highest doubles of the real parts of Q^H*b;
 *   brelohihi are the second highest doubles of the real parts of Q^H*b;
 *   brehilohi are the third highest doubles of the real parts of Q^H*b;
 *   brelolohi are the fourth highest doubles of the real parts of Q^H*b;
 *   brehihilo are the fourth lowest doubles of the real parts of Q^H*b;
 *   brelohilo are the third lowest doubles of the real parts of Q^H*b;
 *   brehilolo are the second lowest doubles of the real parts of Q^H*b;
 *   brelololo are the lowest doubles of the real parts of Q^H*b;
 *   bimhihihi are the highest doubles of the imaginary parts of Q^H*b;
 *   bimlohihi are the second highest doubles of the imaginary parts of Q^H*b;
 *   bimhilohi are the third highest doubles of the imaginary parts of Q^H*b;
 *   bimlolohi are the fourth highest doubles of the imaginary parts of Q^H*b;
 *   bimhihilo are the fourth lowest doubles of the imaginary parts of Q^H*b;
 *   bimlohilo are the third lowest doubles of the imaginary parts of Q^H*b;
 *   bimhilolo are the second lowest doubles of the imaginary parts of Q^H*b;
 *   bimlololo are the lowest doubles of the imaginary parts of Q^H*b;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs. */

void GPU_dbl8_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format,
 *   using QR factorization and substitutions, on real data.
 *
 * REQUIRED : dim = szt*nbt.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   tailidx  the index of the start of the update in the tail;
 *   mathihihi are degp1 matrices of dimension dim;
 *   matlohihi are degp1 matrices of dimension dim;
 *   mathilohi are degp1 matrices of dimension dim;
 *   matlolohi are degp1 matrices of dimension dim;
 *   mathihilo are degp1 matrices of dimension dim;
 *   matlohilo are degp1 matrices of dimension dim;
 *   mathilolo are degp1 matrices of dimension dim;
 *   matlololo are degp1 matrices of dimension dim;
 *   Qhihihi  space for the highest doubles of the Q of the QR;
 *   Qlohihi  space for the second highest doubles of the Q of the QR;
 *   Qhilohi  space for the third highest doubles of the Q of the QR;
 *   Qlolohi  space for the fourth highest doubles of the Q of the QR;
 *   Qhihilo  space for the fourth lowest doubles of the Q of the QR;
 *   Qlohilo  space for the third lowest doubles of the Q of the QR;
 *   Qhilolo  space for the second lowest doubles of the Q of the QR;
 *   Qlololo  space for the lowest doubles of the Q of the QR;
 *   Rhihihi  space for the highest doubles of the R of the QR;
 *   Rlohihi  space for the second highest doubles of the R of the QR;
 *   Rhilohi  space for the third highest doubles of the R of the QR;
 *   Rlolohi  space for the fourth highest doubles of the R of the QR;
 *   Rhihilo  space for the fourth lowst doubles of the R of the QR;
 *   Rlohilo  space for the third lowest doubles of the R of the QR;
 *   Rhilolo  space for the second lowest doubles of the R of the QR;
 *   Rlololo  space for the lowest doubles of the R of the QR;
 *   rhshihihi are degp1 vectors of dimension dim;
 *   rhslohihi are degp1 vectors of dimension dim;
 *   rhshilohi are degp1 vectors of dimension dim;
 *   rhslolohi are degp1 vectors of dimension dim;
 *   rhshihilo are degp1 vectors of dimension dim;
 *   rhslohilo are degp1 vectors of dimension dim;
 *   rhshilolo are degp1 vectors of dimension dim;
 *   rhslololo are degp1 vectors of dimension dim;
 *   solhihihi has space allocated for degp1 vectors of dimension dim;
 *   sollohihi has space allocated for degp1 vectors of dimension dim;
 *   solhilohi has space allocated for degp1 vectors of dimension dim;
 *   sollolohi has space allocated for degp1 vectors of dimension dim;
 *   solhihilo has space allocated for degp1 vectors of dimension dim;
 *   sollohilo has space allocated for degp1 vectors of dimension dim;
 *   solhilolo has space allocated for degp1 vectors of dimension dim;
 *   sollololo has space allocated for degp1 vectors of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qhihihi  highest doubles of the Q of the QR of the Jacobian;
 *   Qlohihi  second highest doubles of the Q of the QR of the Jacobian;
 *   Qhilohi  third highest doubles of the Q of the QR of the Jacobian;
 *   Qlolohi  fourth highest doubles of the Q of the QR of the Jacobian;
 *   Qhihilo  fourth lowest doubles of the Q of the QR of the Jacobian;
 *   Qlohilo  third lowest doubles of the Q of the QR of the Jacobian;
 *   Qhilolo  second lowest doubles of the Q of the QR of the Jacobian;
 *   Qlololo  lowest doubles of the Q of the QR of the Jacobian;
 *   Rhihihi  highest doubles of the R of the QR of the Jacobian;
 *   Rlohihi  second highest doubles of the R of the QR of the Jacobian;
 *   Rhilohi  third highest doubles of the R of the QR of the Jacobian;
 *   Rlolohi  fourth highest doubles of the R of the QR of the Jacobian;
 *   Rhihilo  fourth lowest doubles of the R of the QR of the Jacobian;
 *   Rlohilo  third lowest doubles of the R of the QR of the Jacobian;
 *   Rhilolo  second lowest doubles of the R of the QR of the Jacobian;
 *   Rlololo  lowest doubles of the R of the QR of the Jacobian;
 *   rhshihihi are the highest doubles of updated right hand sides;
 *   rhslohihi are the second highest doubles of updated right hand sides;
 *   rhshilohi are the third highest doubles of updated right hand sides;
 *   rhslolohi are the fourth highest doubles of updated right hand sides;
 *   rhshihilo are the fourth lowest doubles of updated right hand sides;
 *   rhslohilo are the third lowest doubles of updated right hand sides;
 *   rhshilolo are the second lowest doubles of updated right hand sides;
 *   rhslololo are the lowest doubles of updated right hand sides;
 *   solhihihi are the highest doubles of the solution series;
 *   sollohihi are the second highest doubles of the solution series;
 *   solhilohi are the third highest doubles of the solution series;
 *   sollolohi are the fourth highest doubles of the solution series;
 *   solhihilo are the fourth lowest doubles of the solution series;
 *   sollohilo are the third lowest doubles of the solution series;
 *   solhilolo are the second lowest doubles of the solution series;
 *   sollololo are the lowest doubles of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates. */

void GPU_cmplx8_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi, 
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo, 
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi, 
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo, 
   double **solimhilolo, double **solimlololo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsed, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format,
 *   using QR factorization and substitutions, on complex data.
 *
 * REQUIRED : dim = szt*nbt.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   szt      size of each block (and tile);
 *   nbt      number of blocks (and tiles) dim = szt*nbt; 
 *   tailidx  the index of the start of the update in the tail;
 *   matrehihihi are degp1 matrices of dimension dim;
 *   matrelohihi are degp1 matrices of dimension dim;
 *   matrehilohi are degp1 matrices of dimension dim;
 *   matrelolohi are degp1 matrices of dimension dim;
 *   matrehihilo are degp1 matrices of dimension dim;
 *   matrelohilo are degp1 matrices of dimension dim;
 *   matrehilolo are degp1 matrices of dimension dim;
 *   matrelololo are degp1 matrices of dimension dim;
 *   matimhihihi are degp1 matrices of dimension dim;
 *   matimlohihi are degp1 matrices of dimension dim;
 *   matimhilohi are degp1 matrices of dimension dim;
 *   matimlolohi are degp1 matrices of dimension dim;
 *   matimhihilo are degp1 matrices of dimension dim;
 *   matimlohilo are degp1 matrices of dimension dim;
 *   matimhilolo are degp1 matrices of dimension dim;
 *   matimlololo are degp1 matrices of dimension dim;
 *   Qrehihihi has space for the highest doubles of
 *            the real parts of the Q of the QR;
 *   Qrelohihi has space for the second highest doubles of
 *            the real parts of the Q of the QR;
 *   Qrehilohi has space for the third highest doubles of
 *            the real parts of the Q of the QR;
 *   Qrelolohi has space for the fourth highest doubles of
 *            the real parts of the Q of the QR;
 *   Qrehihilo has space for the fourth lowest doubles of 
 *            the real parts of the Q of the QR;
 *   Qrelohilo has space for the third lowest doubles of 
 *            the real parts of the Q of the QR;
 *   Qrehilolo has space for the second lowest doubles of 
 *            the real parts of the Q of the QR;
 *   Qrelololo has space for the lowest doubles of 
 *            the real parts of the Q of the QR;
 *   Qimhihihi has  space for the highest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimlohihi has  space for the second highest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimhilohi has space for the third highest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimlolohi has space for the fourth highest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimhihilo has space for the fourth lowest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimlohilo has space for the third lowest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimhilolo has space for the second lowest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Qimlololo has space for the lowest doubles of
 *            the imaginary parts of the Q of the QR;
 *   Rrehihihi has space for the highest doubles of
 *            the real parts of the R of the QR;
 *   Rrelohihi has space for the second highest doubles of
 *            the real parts of the R of the QR;
 *   Rrehilohi has space for the third highest doubles of
 *            the real parts of the R of the QR;
 *   Rrelolohi has space for the fourth highest doubles of
 *            the real parts of the R of the QR;
 *   Rrehihilo has space for the fourth lowest doubles of 
 *            the real parts of the R of the QR;
 *   Rrelohilo has space for the third lowest doubles of 
 *            the real parts of the R of the QR;
 *   Rrehilolo has space for the second lowest doubles of 
 *            the real parts of the R of the QR;
 *   Rrelololo has space for the lowest doubles of 
 *            the real parts of the R of the QR;
 *   Rimhihihi has space for the highest doubles of
 *            the imaginary parts of the R of the QR;
 *   Rimlohihi has space for the second highest doubles of
 *            the imaginary parts of the R of the QR;
 *   Rimhilohi has space for the third highest doubles of
 *            the imaginary parts of the R of the QR;
 *   Rimlolohi has space for the fourth highest doubles of
 *            the imaginary parts of the R of the QR;
 *   Rimhihilo has space for the fourth lowest doubles of 
 *            the imaginary parts of the R of the QR;
 *   Rimlohilo has space for the third lowest doubles of 
 *            the imaginary parts of the R of the QR;
 *   Rimhilolo has space for the second lowest doubles of 
 *            the imaginary parts of the R of the QR;
 *   Rimlololo has space for the lowest doubles of 
 *            the imaginary parts of the R of the QR;
 *   rhsrehihihi are degp1 vectors of dimension dim;
 *   rhsrelohihi are degp1 vectors of dimension dim;
 *   rhsrehilohi are degp1 vectors of dimension dim;
 *   rhsrelolohi are degp1 vectors of dimension dim;
 *   rhsrehihilo are degp1 vectors of dimension dim;
 *   rhsrelohilo are degp1 vectors of dimension dim;
 *   rhsrehilolo are degp1 vectors of dimension dim;
 *   rhsrelololo are degp1 vectors of dimension dim;
 *   rhsimhihihi are degp1 vectors of dimension dim;
 *   rhsimlohihi are degp1 vectors of dimension dim;
 *   rhsimhilohi are degp1 vectors of dimension dim;
 *   rhsimlolohi are degp1 vectors of dimension dim;
 *   rhsimhihilo are degp1 vectors of dimension dim;
 *   rhsimlohilo are degp1 vectors of dimension dim;
 *   rhsimhilolo are degp1 vectors of dimension dim;
 *   rhsimlololo are degp1 vectors of dimension dim;
 *   solrehihihi has space allocated for degp1 vectors of dimension dim;
 *   solrelohihi has space allocated for degp1 vectors of dimension dim;
 *   solrehilohi has space allocated for degp1 vectors of dimension dim;
 *   solrelolohi has space allocated for degp1 vectors of dimension dim;
 *   solrehihilo has space allocated for degp1 vectors of dimension dim;
 *   solrelohilo has space allocated for degp1 vectors of dimension dim;
 *   solrehilolo has space allocated for degp1 vectors of dimension dim;
 *   solrelololo has space allocated for degp1 vectors of dimension dim;
 *   solimhihihi has space allocated for degp1 vectors of dimension dim;
 *   solimlohihi has space allocated for degp1 vectors of dimension dim;
 *   solimhilohi has space allocated for degp1 vectors of dimension dim;
 *   solimlolohi has space allocated for degp1 vectors of dimension dim;
 *   solimhihilo has space allocated for degp1 vectors of dimension dim;
 *   solimlohilo has space allocated for degp1 vectors of dimension dim;
 *   solimhilolo has space allocated for degp1 vectors of dimension dim;
 *   solimlololo has space allocated for degp1 vectors of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   Qrehihihi are the highest doubles of the real parts of Q;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third highest doubles of the real parts of Q;
 *   Qrelolohi are the fourth highest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo are the third lowest doubles of the real parts of Q;
 *   Qrehilolo are the second lowest doubles of the real parts of Q;
 *   Qrelololo are the lowest doubles of the real parts of Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi are the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q;
 *   Rrehihihi are the highest doubles of the real parts of R;
 *   Rrelohihi are the second highest doubles of the real parts of R;
 *   Rrehilohi are the third highest doubles of the real parts of R;
 *   Rrelolohi are the fourth highest doubles of the real parts of R;
 *   Rrehihilo are the fourth lowest doubles of the real parts of R;
 *   Rrelohilo are the third lowest doubles of the real parts of R;
 *   Rrehilolo are the second lowest doubles of the real parts of R;
 *   Rrelololo are the lowest doubles of the real parts of R;
 *   Rimhihihi are the highest doubles of the imaginary parts of R;
 *   Rimlohihi are the second highest doubles of the imaginary parts of R;
 *   Rimhilohi are the third highest doubles of the imaginary parts of R;
 *   Rimlolohi are the fourth highest doubles of the imaginary parts of R;
 *   Rimhihilo are the fourth lowest doubles of the imaginary parts of R;
 *   Rimlohilo are the third lowest doubles of the imaginary parts of R;
 *   Rimhilolo are the second lowest doubles of the imaginary parts of R;
 *   Rimlololo are the lowest doubles of the imaginary parts of R;
 *   rhsrehihi are the highest doubles of the real parts of rhs;
 *   rhsrelohihi are the 2nd highest doubles of the real parts of rhs;
 *   rhsrehilohi are the 3rd highest doubles of the real parts of rhs;
 *   rhsrelolohi are the 4th highest doubles of the real parts of rhs;
 *   rhsrehihilo are the 4th lowest doubles of the real parts of rhs;
 *   rhsrelohilo are the 3rd lowest doubles of the real parts of rhs;
 *   rhsrehilolo are the 2nd lowest doubles of the real parts of rhs;
 *   rhsrelololo are the lowest doubles of the real parts of rhs;
 *   rhsimhihihi are the highest doubles of the imag parts of rhs;
 *   rhsimlohihi are the 2nd highest doubles of the imag parts of rhs;
 *   rhsimhilohi are the 3rd highest doubles of the imag parts of rhs;
 *   rhsimlolohi are the 4th highest doubles of the imag parts of rhs;
 *   rhsimhihilo are the 4th lowest doubles of the imag parts of rhs;
 *   rhsimlohilo are the 3rd lowest doubles of the imag parts of rhs;
 *   rhsimhilolo are the 2nd lowest doubles of the imag parts of rhs;
 *   rhsimlololo are the lowest doubles of the imag parts of rhs;
 *   solrehihihi are the 2nd highest doubles of the real parts of the solution;
 *   solrelohihi are the 3rd highest doubles of the real parts of the solution;
 *   solrehilohi are the 4th highest doubles of the real parts of the solution;
 *   solrehihilo are the 4th lowest doubles of the real parts of the solution;
 *   solrelohilo are the 3rd lowest doubles of the real parts of the solution;
 *   solrehilolo are the 2nd lowest doubles of the real parts of the solution;
 *   solrelololo are the lowest doubles of the real parts of the solution;
 *   solimhihihi are the highest doubles of the imag parts of the solution;
 *   solimlohihi are the 2nd highest doubles of the imag parts of the solution;
 *   solimhilohi are the 3rd highest doubles of the imag parts of the solution;
 *   solimlolohi are the 4th highest doubles of the imag parts of the solution;
 *   solimhihilo are the 4th lowest doubles of the imag parts of the solution;
 *   solimlohilo are the 3rd lowest doubles of the imag parts of the solution;
 *   solimhilolo are the 2nd lowest doubles of the imag parts of the solution;
 *   solimlololo are the lowest doubles of the imag parts of the solution;
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
