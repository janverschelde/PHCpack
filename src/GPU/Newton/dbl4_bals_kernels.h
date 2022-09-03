// The file dbl4_bals_kernels.h defines the prototypes of functions
// to define memory transfers and kernel launches to solve a linear system
// of power series linearization, in quad double precision.

#ifndef __dbl4_bals_kernels_h__
#define __dbl4_bals_kernels_h__

__global__ void dbl4_bals_tail
 ( int ncols, int szt,
   double *Ahihi, double *Alohi, double *Ahilo, double *Alolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo );
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
 *   Ahihi    nrows-by-ncols matrix of highest doubles to multiply x with;
 *   Alohi    nrows-by-ncols matrix of second highest doubles;
 *   Ahilo    nrows-by-ncols matrix of second lowest doubles;
 *   Alolo    nrows-by-ncols matrix of lowest doubles;
 *   xhihi    highest doubles of vector of dimension ncols;
 *   xlohi    second highest doubles of vector of dimension ncols;
 *   xhilo    second lowest doubles of vector of dimension ncols;
 *   xlolo    lowest doubles of vector of dimension ncols;
 *   bhihi    highest doubles of vector of dimension nrows;
 *   blohi    second highest doubles of vector of dimension nrows.
 *   bhilo    second lowest doubles of vector of dimension nrows;
 *   blolo    lowest doubles of vector of dimension nrows.
 *
 * ON RETURN :
 *   bhihi    highest doubles of the updated right hand side;
 *   blohi    second highest doubles of the updated right hand side;
 *   bhilo    second lowest doubles of the updated right hand side;
 *   blolo    lowest doubles of the updated right hand side. */

__global__ void dbl4_bals_qtb
 ( int ncols, int szt,
   double *Qthihi, double *Qtlohi, double *Qthilo, double *Qtlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *rhihi, double *rlohi, double *rhilo, double *rlolo );
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

void GPU_dbl4_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo, 
   double *xhihi, double *xlohi, double *xhilo, double *xlolo, bool verbose );
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
 *   verbose  is the verbose flag.
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
 *   xlolo    lowest doubles of the least squares solution. */

void GPU_dbl4_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
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

void GPU_dbl4_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo, bool verbose );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   bhihi    highest doubles of the product of Q^T with b;
 *   blohi    second highest doubles of the product of Q^T with b;
 *   bhilo    second lowest doubles of the product of Q^T with b;
 *   blolo    lowest doubles of the product of Q^T with b. */

void GPU_dbl4_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
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
 *   vrblvl   the verbose level (0 for silent).
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
 *   sollolo  lowest doubles of the solution series. */

#endif
