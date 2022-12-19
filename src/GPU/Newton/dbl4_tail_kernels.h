// The file dbl4_tail_kernels.h defines the prototypes of functions to
// define memory transfers and kernel launches for the tail operations to
// solve a linear system of power series, in quad double precision.

#ifndef __dbl4_tail_kernels_h__
#define __dbl4_tail_kernels_h__

__global__ void dbl4_bals_tail
 ( int ncols, int szt,
   double *Ahihi, double *Alohi, double *Ahilo, double *Alolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo );
/*
 * DESCRIPTION :
 *   Subtracts from the right hand side b the product of A with x,
 *   on real data.
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

__global__ void cmplx4_bals_tail
 ( int ncols, int szt,
   double *Arehihi, double *Arelohi, double *Arehilo, double *Arelolo,
   double *Aimhihi, double *Aimlohi, double *Aimhilo, double *Aimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo, 
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo );
/*
 * DESCRIPTION :
 *   Subtracts from the right hand side b the product of A with x,
 *   on complex data.
 *
 * REQUIRED : nrows = szt times the number of blocks,
 *   where nrows in the number of rows in A and the dimension of b.
 *
 * ON ENTRY :
 *   ncols    number of columns in A and the dimension of x;
 *   szt      size of each block (and tile);
 *   Arehihi  highest doubles of the real parts of nrows-by-ncols matrix;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the real parts of A;
 *   Arelolo  lowest doubles of the real parts of A;
 *   Aimhihi  highest doubles of the imaginary parts of A;
 *   Aimlohi  second highest doubles of the imaginary parts of A;
 *   Aimhilo  second lowest doubles of the imaginary parts of A;
 *   Aimlolo  lowest doubles of the imaginary parts of A;
 *   xrehihi  highest doubles of the real parts of vector of dimension ncols;
 *   xrelohi  second highest doubles of the real parts of x;
 *   xrehilo  second lowest doubles of the real parts of x;
 *   xrelolo  lowest doubles of the real parts of x;
 *   ximhihi  highest doubles of the imaginary parts of x; 
 *   ximlohi  second highest doubles of the imaginary parts of x; 
 *   ximhilo  second lowest doubles of the imaginary parts of x;
 *   ximlolo  lowest doubles of the imaginary parts of x;
 *   brehihi  highest doubles of the real parts of vector of dimension nrows;
 *   brelohi  second highest doubles of the real parts of b;
 *   brehilo  second lowest doubles of the real parts of b;
 *   brelolo  lowest doubles of the real parts of b;
 *   bimhihi  highest doubles of the imaginary parts of b;
 *   bimlohi  second highest doubles of the imaginary parts of b;
 *   bimhilo  second lowest doubles of the imaginary parts of b;
 *   bimlolo  lowest doubles of the imaginary parts of b.
 *
 * ON RETURN :
 *   brehihi  highest doubles of the real parts of updated right hand side;
 *   brelohi  second highest doubles of the real parts of b;
 *   brehilo  second lowest doubles of the real parts of b;
 *   brelolo  lowest doubles of the real parts of b;
 *   bimhihi  highest doubles of the imag parts of b;
 *   bimlohi  second highest doubles of the imag parts of b;
 *   bimhilo  second lowest doubles of the imag parts of b;
 *   bimlolo  lowest doubles of the imag parts of b. */

void write_dbl4_balsflops ( int ctype, int ncols, float lapsms );
/*
 * DESCRIPTION :
 *   Writes the number of flops and arithmetic intensity for b = b - A*x.
 *
 * ON ENTRY :
 *   ctype    0 if real, 1 if complex;
 *   ncols    number of columns in the matrix;
 *   lapsms   time elapsed in milliseconds. */ 

void GPU_dbl4_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double *totupdlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   After each block of coefficients of the series,
 *   kernels are launched for the multiplication of the tail matrices
 *   with the solution coefficients to update the right hand sides
 *   of the linear system of power series, on real data.
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
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   rhshihi  highest doubles of updated right hand sides;
 *   rhslohi  second highest doubles of updated right hand sides;
 *   rhshilo  second lowest doubles of updated right hand sides;
 *   rhslolo  lowest doubles of updated right hand sides;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds. */

void GPU_cmplx4_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double *totupdlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   After each block of coefficients of the series,
 *   kernels are launched for the multiplication of the tail matrices
 *   with the solution coefficients to update the right hand sides
 *   of the linear system of power series, on complex data.
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
 *   matrehihi are the highest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrelohi are the second highest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrehilo are the second lowest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrelolo are the lowest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matimhihi are the highest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimlohi are the second highest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimhilo are the second lowest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimlolo are the lowest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   rhsrehihi are the highest doubles of the real parts of rhs;
 *   rhsrelohi are the second highest doubles of the real parts of rhs;
 *   rhsrehilo are the second lowest doubles of the real parts of rhs;
 *   rhsrelolo are the lowest doubles of the real parts of rhs;
 *   rhsimhihi are the highest doubles of the imaginary parts of rhs;
 *   rhsimlohi are the second highest doubles of the imaginary parts of rhs;
 *   rhsimhilo are the second lowest doubles of the imaginary parts of rhs;
 *   rhsimlolo are the lowest doubles of the imaginary parts of rhs;
 *   solrehihi are the highest doubles of the real parts of the solution,
 *            computed up to stage-1;
 *   solrelohi are the second highest doubles of the real parts of the
 *            solution, computed up to stage-1;
 *   solrehilo are the second lowest doubles of the real parts of the
 *            solution, computed up to stage-1;
 *   solrelolo are the lowest doubles of the real parts of the solution,
 *            computed up to stage-1;
 *   solimhihi are the highest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimlohi are the second highest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimhilo are the second lowest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimlolo are the lowest doubles of the imaginary parts of the solution,
 *            computed up to stage-1;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   rhsrehihi are the highest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrelohi are the second highest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrehilo are the second lowest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrelolo are the lowest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsimhihi are the highest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimlohi are the second highest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimhilo are the second lowest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimlolo are the lowest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds. */

void GPU_dbl4_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **resvechihi, double **resveclohi,
   double **resvechilo, double **resveclolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
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
 *   mathihi  degp1 matrices of dimension dim;
 *   matlohi  degp1 matrices of dimension dim;
 *   mathilo  degp1 matrices of dimension dim;
 *   matlolo  degp1 matrices of dimension dim;
 *   rhshihi  degp1 right hand side vectors of dimension dim;
 *   rhslohi  degp1 right hand side vectors of dimension dim;
 *   rhshilo  degp1 right hand side vectors of dimension dim;
 *   rhslolo  degp1 right hand side vectors of dimension dim;
 *   solhihi  degp1 solution vectors of dimension dim;
 *   sollohi  degp1 solution vectors of dimension dim;
 *   solhilo  degp1 solution vectors of dimension dim;
 *   sollolo  degp1 solution vectors of dimension dim;
 *   resvechihi has space for the residual power series;
 *   resveclohi has space for the residual power series;
 *   resvechilo has space for the residual power series;
 *   resveclolo has space for the residual power series;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   resvechihi are the highest doubles of the residual power series;
 *   resveclohi are the 2nd highest doubles of the residual power series;
 *   resvechilo are the 2nd lowest doubles of the residual power series;
 *   resveclolo are the lowest doubles of the residual power series;
 *   resmaxhihi is the highest double of the maximum component in resvec;
 *   resmaxlohi is the 2nd highest double of the maximum component in resvec;
 *   resmaxhilo is the 2nd lowest double of the maximum component in resvec;
 *   resmaxlolo is the lowest double of the maximum component in resvec;
 *   lapms    elapsed time spent by all kernels, in milliseconds;
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx4_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi, 
   double **rhsimhilo, double **rhsimlolo, 
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
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
 *   matrehihi are degp1 matrices of dimension dim;
 *   matrelohi are degp1 matrices of dimension dim;
 *   matrehilo are degp1 matrices of dimension dim;
 *   matrelolo are degp1 matrices of dimension dim;
 *   matimhihi are degp1 matrices of dimension dim;
 *   matimlohi are degp1 matrices of dimension dim;
 *   matimhilo are degp1 matrices of dimension dim;
 *   matimlolo are degp1 matrices of dimension dim;
 *   rhsrehihi are degp1 right hand side vectors of dimension dim;
 *   rhsrelohi are degp1 right hand side vectors of dimension dim;
 *   rhsrehilo are degp1 right hand side vectors of dimension dim;
 *   rhsrelolo are degp1 right hand side vectors of dimension dim;
 *   rhsimhihi are degp1 right hand side vectors of dimension dim;
 *   rhsimlohi are degp1 right hand side vectors of dimension dim;
 *   rhsimhilo are degp1 right hand side vectors of dimension dim;
 *   rhsimlolo are degp1 right hand side vectors of dimension dim;
 *   solrehihi are degp1 solution vectors of dimension dim;
 *   solrelohi are degp1 solution vectors of dimension dim;
 *   solrehilo are degp1 solution vectors of dimension dim;
 *   solrelolo are degp1 solution vectors of dimension dim;
 *   solimhihi are degp1 solution vectors of dimension dim;
 *   solimlohi are degp1 solution vectors of dimension dim;
 *   solimhilo are degp1 solution vectors of dimension dim;
 *   solimlolo are degp1 solution vectors of dimension dim;
 *   resvecrehihi has space for the residual power series;
 *   resvecrelohi has space for the residual power series;
 *   resvecrehilo has space for the residual power series;
 *   resvecrelolo has space for the residual power series;
 *   resvecimhihi has space for the residual power series;
 *   resvecimlohi has space for the residual power series;
 *   resvecimhilo has space for the residual power series;
 *   resvecimlolo has space for the residual power series;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   resvecrehihi are the highest doubles of the real parts
 *            of the residual series;
 *   resvecrelohi are the second highest doubles of the real parts
 *            of the residual series;
 *   resvecrehilo are the second lowest doubles of the real parts
 *            of the residual series;
 *   resvecrelolo are the lowest doubles of the real parts
 *            of the residual series;
 *   resvecimhihi are the highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlohi are the second highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimhilo are the second lowest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlolo are the lowest doubles of the imaginary parts
 *            of the residual series;
 *   resmaxhihi is the highest double of the max norm of the residual;
 *   resmaxlohi is the second highest double of the max norm of the residual;
 *   resmaxhilo is the second lowest double of the max norm of the residual;
 *   resmaxlolo is the lowest double of the max norm of the residual;
 *   lapms    elapsed time spent by all kernels, in milliseconds;
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

#endif
