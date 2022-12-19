// The file dbl_tail_kernels.h defines the prototypes of functions to
// define memory transfers and kernel launches for the tail operations to
// solve a linear system of power series, in double precision.

#ifndef __dbl_tail_kernels_h__
#define __dbl_tail_kernels_h__

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

void write_dbl_balsflops ( int ctype, int ncols, float lapsms );
/*
 * DESCRIPTION :
 *   Writes the number of flops and arithmetic intensity for b = b - A*x.
 *
 * ON ENTRY :
 *   ctype    0 if real, 1 if complex;
 *   ncols    number of columns in the matrix;
 *   lapsms   time elapsed in milliseconds. */ 

void GPU_dbl_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mat, double **rhs, double **sol,
   double *totupdlapsedms, int vrblvl );
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
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   rhs      updated right hand sides;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds. */

void GPU_cmplx_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matre, double ***matim, double **rhsre, double **rhsim,
   double **solre, double **solim, double *totupdlapsedms, int vrblvl );
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
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   rhsre    real parts of the updated right hand sides;
 *   rhsim    imaginary parts of the updated right hand sides;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds. */

void GPU_dbl_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mat, double **rhs, double **sol, double **resvec, double *resmax,
   double *totreslapsedms, long long int *add, long long int *mul,
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
 *   tailidx  the index of the start of the update in the tail;
 *   mat      degp1 matrices of dimension dim;
 *   rhs      degp1 right hand side vectors of dimension dim;
 *   sol      degp1 solution vectors of dimension dim;
 *   resvec   space for the residual power series;
 *   totreslapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   resvec   the residual power series;
 *   resmax   maximum component of the residual power series;
 *   totreslapsedms acculumates time spent by all kernels in milliseconds;
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matre, double ***matim, double **rhsre, double **rhsim,
   double **solre, double **solim,
   double **resvecre, double **resvecim, double *resmax,
   double *totreslapsedms, long long int *add, long long int *mul,
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
 *   tailidx  the index of the start of the update in the tail;
 *   matre    degp1 matrices of dimension dim;
 *   matim    degp1 matrices of dimension dim;
 *   rhsre    degp1 right hand side vectors of dimension dim;
 *   rhsim    degp1 right hand side vectors of dimension dim;
 *   solre    degp1 solution vectors of dimension dim;
 *   solim    degp1 solution vectors of dimension dim;
 *   resvecre has space for the residual power series;
 *   resvecim has space for the residual power series;
 *   totreslapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   resvecre are the real parts of the residual power series;
 *   resvecim are the imaginary parts the residual power series;
 *   resmax   max norm of the residual power series;
 *   totreslapsedms acculumates time spent by all kernels in milliseconds;
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

#endif
