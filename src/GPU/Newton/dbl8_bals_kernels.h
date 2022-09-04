// The file dbl8_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in octo double precision.

#ifndef __dbl8_bals_kernels_h__
#define __dbl8_bals_kernels_h__

__global__ void dbl8_bals_tail
 ( int ncols, int szt,
   double *Ahihihi, double *Alohihi, double *Ahilohi, double *Alolohi,
   double *Ahihilo, double *Alohilo, double *Ahilolo, double *Alololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo );
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
 *   Ahihihi  nrows-by-ncols matrix of highest doubles of A;
 *   Alohihi  rows-by-ncols matrix of second highest doubles;
 *   Ahilohi  nrows-by-ncols matrix of second highest doubles;
 *   Alolohi  nrows-by-ncols matrix of fourth highest doubles;
 *   Ahihilo  nrows-by-ncols matrix of fourth lowest doubles of A;
 *   Alohilo  nrows-by-ncols matrix of third lowest doubles;
 *   Ahilolo  nrows-by-ncols matrix of second lowest doubles;
 *   Alololo  nrows-by-ncols matrix of lowest doubles;
 *   xhihihi  highest doubles of vector of dimension ncols;
 *   xlohihi  second highest doubles of vector of dimension ncols;
 *   xhilohi  third highest doubles of vector of dimension ncols;
 *   xlolohi  fourth highest doubles of vector of dimension ncols;
 *   xhihilo  fourth lowest doubles of vector of dimension ncols;
 *   xlohilo  third lowest doubles of vector of dimension ncols;
 *   xhilolo  second lowest doubles of vector of dimension ncols;
 *   xlololo  lowest doubles of vector of dimension ncols;
 *   bhihihi  highest doubles of vector of dimension nrows;
 *   blohihi  second highest doubles of vector of dimension nrows;
 *   bhilohi  third highest doubles of vector of dimension nrows;
 *   blolohi  fourth highest doubles of vector of dimension nrows;
 *   bhihilo  fourth lowest doubles of vector of dimension nrows;
 *   blohilo  third lowest doubles of vector of dimension nrows;
 *   bhilolo  second lowest doubles of vector of dimension nrows;
 *   blololo  lowest doubles of vector of dimension nrows.
 *
 * ON RETURN :
 *   bhihihi  highest doubles of the updated right hand side;
 *   blohihi  second highest doubles of the updated right hand side;
 *   bhilohi  third highest doubles of the updated right hand side;
 *   blolohi  fourth highest doubles of the updated right hand side;
 *   bhihilo  fourth lowest doubles of the updated right hand side;
 *   blohilo  third lowest doubles of the updated right hand side;
 *   bhilolo  second lowest doubles of the updated right hand side;
 *   blololo  lowest doubles of the updated right hand side. */

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
 *   Multiplies the transpose of Q with b.
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
   bool verbose );
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
 *   verbose  is the verbose flag.
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
 *   xlololo  lowest doubles of the least squares solution. */

void GPU_dbl8_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo, bool verbose );
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
 *   mathihi  highest doubles of the matrices of the linearized series;
 *   matlohi  second highest doubles of the matrices of the linearized series;
 *   mathilo  second lowest doubles of the matrices of the linearized series;
 *   matlolo  lowest doubles of the matrices of the linearized series;
 *   rhshihi  highest doubles of the right hand sides of the system;
 *   rhslohi  second highest doubles of the right hand sides of the system;
 *   rhshilo  second lowest doubles of the right hand sides of the system;
 *   rhslolo  lowest doubles of the right hand sides of the system;
 *   solhihi  highest doubles of solution computed up to stage-1;
 *   sollohi  second highest doubles of solution computed up to stage-1;
 *   solhilo  second lowest doubles of solution computed up to stage-1;
 *   sollolo  lowest doubles of solution computed up to stage-1;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   rhshihi  highest doubles of updated right hand sides;
 *   rhslohi  second highest doubles of updated right hand sides;
 *   rhshilo  second lowest doubles of updated right hand sides;
 *   rhslolo  lowest doubles of updated right hand sides. */

void GPU_dbl8_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   bool verbose );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   bhihihi  highest doubles of the product of Q^T with b;
 *   blohihi  second highest doubles of the product of Q^T with b;
 *   bhilohi  third highest doubles of the product of Q^T with b;
 *   blolohi  fourth highest doubles of the product of Q^T with b;
 *   bhihilo  fourth lowest doubles of the product of Q^T with b;
 *   blohilo  third lowest doubles of the product of Q^T with b;
 *   bhilolo  second lowest doubles of the product of Q^T with b;
 *   blololo  lowest doubles of the product of Q^T with b. */

void GPU_dbl8_bals_solve
 ( int dim, int degp1, int szt, int nbt,
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
   double **solhilolo, double **sollololo, int vrblvl );
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
 *   vrblvl   the verbose level (0 for silent).
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
 *   sollololo are the lowest doubles of the solution series. */

#endif
