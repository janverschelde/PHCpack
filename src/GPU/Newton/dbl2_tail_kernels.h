// The file dbl2_tail_kernels.h defines the prototypes of functions to
// define memory transfers and kernel launches for the tail operations ot
// solve a linear system of power series, in double double precision.

#ifndef __dbl2_tail_kernels_h__
#define __dbl2_tail_kernels_h__

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

void write_dbl2_balsflops ( int ctype, int ncols, float lapsms );
/*
 * DESCRIPTION :
 *   Writes the number of flops and arithmetic intensity for b = b - A*x.
 *
 * ON ENTRY :
 *   ctype    0 if real, 1 if complex;
 *   ncols    number of columns in the matrix;
 *   lapsms   time elapsed in milliseconds. */ 

void GPU_dbl2_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathi, double ***matlo, double **rhshi, double **rhslo,
   double **solhi, double **sollo, double *totupdlapsedms, int vrblvl );
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
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   rhshi    high doubles of updated right hand sides;
 *   rhslo    low doubles of updated right hand sides;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds. */

void GPU_cmplx2_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
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
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   rhsrehi  high doubles of the real parts of the updated right hand sides;
 *   rhsrelo  low doubles of the real parts of the updated right hand sides;
 *   rhsimhi  high doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimlo  low doubles of the imaginary parts
 *            of the updated right hand sides;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds. */

void GPU_dbl2_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   double *lapms, long long int *add, long long int *mul, int vrblvl );
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
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 right hand side vectors of dimension dim;
 *   rhslo    degp1 right hand side vectors of dimension dim;
 *   solhi    degp1 solution vectors of dimension dim;
 *   sollo    degp1 solution vectors of dimension dim;
 *   resvechi has space for the residual power series;
 *   resveclo has space for the residual power series;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   resvechi are the high doubles of the residual power series;
 *   resveclo are the low doubles of the residual power series;
 *   resmaxhi is the high double of the maximum component in resvec;
 *   resmaxlo is the low double of the maximum component in resvec;
 *   lapms    elapsed time spent by all kernels, in milliseconds;
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx2_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo, 
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
   double *lapms, long long int *add, long long int *mul, int vrblvl );
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
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   resvecrehi are the high doubles of the real parts of the residual series;
 *   resvecrelo are the low doubles of the real parts of the residual series;
 *   resvecimhi are the high doubles of the imag parts the residual series;
 *   resvecimlo are the low doubles of the imag parts the residual series;
 *   resmaxhi is the high double of the max norm of the residual series;
 *   resmaxlo is the low double of the max norm of the residual series;
 *   lapms    elapsed time spent by all kernels, in milliseconds;
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

#endif
