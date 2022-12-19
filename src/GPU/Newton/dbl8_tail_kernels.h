// The file dbl8_bals_kernels.h defines the prototypes of functions to
// define memory transfers and kernel launches for the tail operations to
// solve a linear system of power series, in octo double precision.

#ifndef __dbl8_tail_kernels_h__
#define __dbl8_tail_kernels_h__

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
 *   Subtracts from the right hand side b the product of A with x,
 *   on real data.
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

__global__ void cmplx8_bals_tail
 ( int ncols, int szt,
   double *Arehihihi, double *Arelohihi, double *Arehilohi, double *Arelolohi,
   double *Arehihilo, double *Arelohilo, double *Arehilolo, double *Arelololo,
   double *Aimhihihi, double *Aimlohihi, double *Aimhilohi, double *Aimlolohi,
   double *Aimhihilo, double *Aimlohilo, double *Aimhilolo, double *Aimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi, 
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo, 
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo );
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
 *   Arehihihi are the highest doubles of the real parts of A;
 *   Arelohihi are the second highest doubles of the real parts of A;
 *   Arehilohi are the third highest doubles of the real parts of A;
 *   Arelolohi are the fourth highest doubles of the real parts of A;
 *   Arehihilo are the fourth lowest doubles of the real parts of A ;
 *   Arelohilo are the third lowest doubles of the real parts of A ;
 *   Arehilolo are the second lowest doubles of the real parts of A ;
 *   Arelololo are the lowest doubles of the real parts of A;
 *   Aimhihihi are the highest doubles of the imaginary parts of A;
 *   Aimlohihi are the second highest doubles of the imaginary parts of A;
 *   Aimhilohi are the third highest doubles of the imaginary parts of A;
 *   Aimlolohi are the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo are the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo are the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo are the second lowest doubles of the imaginary parts of A;
 *   Aimlololo are the lowest doubles of the imaginary parts of A;
 *   xrehihihi are the highest doubles of the real parts of x;
 *   xrelohihi are the second highest doubles of the real parts of x;
 *   xrehilohi are the third highest doubles of the real parts of x;
 *   xrelolohi are the fourth highest doubles of the real parts of x;
 *   xrehihilo are the fourth lowest doubles of the real parts of x;
 *   xrelohilo are the third lowest doubles of the real parts of x;
 *   xrehilolo are the second lowest doubles of the real parts of x;
 *   xrelololo are the lowest doubles of the real parts of x;
 *   ximhihihi are the highest doubles of the imaginary parts of x; 
 *   ximlohihi are the second highest doubles of the imaginary parts of x; 
 *   ximhilohi are the second highest doubles of the imaginary parts of x; 
 *   ximlolohi are the second highest doubles of the imaginary parts of x; 
 *   ximhihilo are the second lowest doubles of the imaginary parts of x;
 *   ximlohilo are the second lowest doubles of the imaginary parts of x;
 *   ximhilolo are the second lowest doubles of the imaginary parts of x;
 *   ximlololo are the lowest doubles of the imaginary parts x;
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
 *   bimlololo are the lowest doubles of the imaginary parts of b.
 *
 * ON RETURN :
 *   brehihihi are the highest doubles of the real parts of updated b;
 *   brelohihi are the second highest doubles of the real parts of b;
 *   brehilohi are the third highest doubles of the real parts of b;
 *   brelolohi are the fourth highest doubles of the real parts of b;
 *   brehihilo are the fourth lowest doubles of the real parts of b;
 *   brelohilo are the third lowest doubles of the real parts of b;
 *   brehilolo are the second lowest doubles of the real parts of b;
 *   brelololo are the lowest doubles of the real parts of b;
 *   bimhihihi are the highest doubles of the imag parts of b;
 *   bimlohihi are the second highest doubles of the imag parts of b;
 *   bimhilohi are the third highest doubles of the imag parts of b;
 *   bimloilohi are the fourth highest doubles of the imag parts of b;
 *   bimhihilo are the fourth lowest doubles of the imag parts of b;
 *   bimlohilo are the third lowest doubles of the imag parts of b;
 *   bimhilolo are the second lowest doubles of the imag parts of b;
 *   bimlololo are the lowest doubles of the imag parts of b. */

void write_dbl8_balsflops ( int ctype, int ncols, float lapsms );
/*
 * DESCRIPTION :
 *   Writes the number of flops and arithmetic intensity for b = b - A*x.
 *
 * ON ENTRY :
 *   ctype    0 if real, 1 if complex;
 *   ncols    number of columns in the matrix;
 *   lapsms   time elapsed in milliseconds. */ 

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
   double **solhilolo, double **sollololo,
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
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds.*/

void GPU_cmplx8_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
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
 *   matrehihihi are the highest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrelohihi are the second highest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrehilohi are the third highest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrelolohi are the fourth highest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrehihilo are the fourth lowest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrelohilo are the third lowest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrehilolo are the second lowest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matrelololo are the lowest doubles of the real parts of
 *            the matrices of the linearized series;
 *   matimhihihi are the highest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimlohihi are the second highest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimhilohi are the third highest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimlolohi are the fourth highest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimhihilo are the fourth lowest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimlohilo are the third lowest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimhilolo are the second lowest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   matimlololo are the lowest doubles of the imaginary parts of
 *            the matrices of the linearized series;
 *   rhsrehihihi are the highest doubles of the real parts of rhs;
 *   rhsrelohihi are the second highest doubles of the real parts of rhs;
 *   rhsrehilohi are the third highest doubles of the real parts of rhs;
 *   rhsrelolohi are the fourth highest doubles of the real parts of rhs;
 *   rhsrehihilo are the fourth lowest doubles of the real parts of rhs;
 *   rhsrelohilo are the third lowest doubles of the real parts of rhs;
 *   rhsrehilolo are the second lowest doubles of the real parts of rhs;
 *   rhsrelololo are the lowest doubles of the real parts of rhs;
 *   rhsimhihihi are the highest doubles of the imaginary parts of rhs;
 *   rhsimlohihi are the second highest doubles of the imaginary parts of rhs;
 *   rhsimhilohi are the third highest doubles of the imaginary parts of rhs;
 *   rhsimlolohi are the fourth highest doubles of the imaginary parts of rhs;
 *   rhsimhihilo are the fourth lowest doubles of the imaginary parts of rhs;
 *   rhsimlohilo are the third lowest doubles of the imaginary parts of rhs;
 *   rhsimhilolo are the second lowest doubles of the imaginary parts of rhs;
 *   rhsimlololo are the lowest doubles of the imaginary parts of rhs;
 *   solrehihihi are the highest doubles of the real parts of the solution,
 *            computed up to stage-1;
 *   solrelohihi are the second highest doubles of the real parts of the
 *            solution, computed up to stage-1;
 *   solrehilohi are the third highest doubles of the real parts of the
 *            solution, computed up to stage-1;
 *   solrelolohi are the fourth highest doubles of the real parts of the
 *            solution, computed up to stage-1;
 *   solrehihilo are the fourth lowest doubles of the real parts of the
 *            solution, computed up to stage-1;
 *   solrelohilo are the third lowest doubles of the real parts of the
 *            solution, computed up to stage-1;
 *   solrehilolo are the second lowest doubles of the real parts of the
 *            solution, computed up to stage-1;
 *   solrelololo are the lowest doubles of the real parts of the solution,
 *            computed up to stage-1;
 *   solimhihihi are the highest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimlohihi are the second highest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimlohihi are the third highest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimhilohi are the fourth highest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimhihilo are the fourth lowest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimlohilo are the third lowest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimhilolo are the second lowest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   solimlololo are the lowest doubles of the imaginary parts of the
 *            solution, computed up to stage-1;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   rhsrehihihi are the highest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrelohihi are the second highest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrehilohi are the third highest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrelolohi are the fourth highest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrehihilo are the fourth lowest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrelohilo are the third lowest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrehilolo are the second lowest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsrelololo are the lowest doubles of the real parts
 *            of the updated right hand sides;
 *   rhsimhihihi are the highest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimlohihi are the second highest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimhilohi are the third highest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimlolohi are the fourth highest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimhihilo are the fourth lowest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimlohilo are the third lowest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimhilolo are the second lowest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   rhsimlololo are the lowest doubles of the imaginary parts
 *            of the updated right hand sides;
 *   totupdlapsedms acculumates time spent by all kernels in milliseconds. */

void GPU_dbl8_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
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
   double **solhilolo, double **sollololo,
   double **resvechihihi, double **resveclohihi,
   double **resvechilohi, double **resveclolohi,
   double **resvechihilo, double **resveclohilo,
   double **resvechilolo, double **resveclololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo,
   double *totreslapsedms, long long int *add, long long int *mul,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear series system, on real data.
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
 *   rhshihihi are degp1 right hand side vectors of dimension dim;
 *   rhslohihi are degp1 right hand side vectors of dimension dim;
 *   rhshilohi are degp1 right hand side vectors of dimension dim;
 *   rhslolohi are degp1 right hand side vectors of dimension dim;
 *   rhshihilo are degp1 right hand side vectors of dimension dim;
 *   rhslohilo are degp1 right hand side vectors of dimension dim;
 *   rhshilolo are degp1 right hand side vectors of dimension dim;
 *   rhslololo are degp1 right hand side vectors of dimension dim;
 *   solhihihi are degp1 solution vectors of dimension dim;
 *   sollohihi are degp1 solution vectors of dimension dim;
 *   solhilohi are degp1 solution vectors of dimension dim;
 *   sollolohi are degp1 solution vectors of dimension dim;
 *   solhihilo are degp1 solution vectors of dimension dim;
 *   sollohilo are degp1 solution vectors of dimension dim;
 *   solhilolo are degp1 solution vectors of dimension dim;
 *   sollololo are degp1 solution vectors of dimension dim;
 *   resvechihihi has space for the residual power series;
 *   resveclohihi has space for the residual power series;
 *   resvechilohi has space for the residual power series;
 *   resveclolohi has space for the residual power series;
 *   resvechihilo has space for the residual power series;
 *   resveclohilo has space for the residual power series;
 *   resvechilolo has space for the residual power series;
 *   resveclololo has space for the residual power series;
 *   totreslapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   resvechihihi are the highest doubles of the residual series;
 *   resveclohihi are the 2nd highest doubles of the residual series;
 *   resvechilohi are the 3rd highest doubles of the residual series;
 *   resveclolohi are the 4th highest doubles of the residual series;
 *   resvechihilo are the 4th lowest doubles of the residual series;
 *   resveclohilo are the 3rd lowest doubles of the residual series;
 *   resvechilolo are the 2nd lowest doubles of the residual series;
 *   resveclololo are the lowest doubles of the residual series;
 *   resmaxhihihi is the highest double of the maximum in resvec;
 *   resmaxlohihi is the 2nd highest double of the maximum in resvec;
 *   resmaxhilohi is the 3rd highest double of the maximum in resvec;
 *   resmaxlolohi is the 4th highest double of the maximum in resvec;
 *   resmaxhihilo is the 4th lowest double of the maximum in resvec;
 *   resmaxlohilo is the 3rd lowest double of the maximum in resvec;
 *   resmaxhilolo is the 2nd lowest double of the maximum in resvec;
 *   resmaxlololo is the lowest double of the maximum in resvec;
 *   totreslapsedms acculumates time spent by all kernels in milliseconds;
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx8_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
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
   double **resvecrehihihi, double **resvecrelohihi,
   double **resvecrehilohi, double **resvecrelolohi,
   double **resvecrehihilo, double **resvecrelohilo,
   double **resvecrehilolo, double **resvecrelololo,
   double **resvecimhihihi, double **resvecimlohihi,
   double **resvecimhilohi, double **resvecimlolohi,
   double **resvecimhihilo, double **resvecimlohilo,
   double **resvecimhilolo, double **resvecimlololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo,
   double *totreslapsedms, long long int *add, long long int *mul,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear series system, on complex data.
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
 *   rhsrehihihi are degp1 right hand side vectors of dimension dim;
 *   rhsrelohihi are degp1 right hand side vectors of dimension dim;
 *   rhsrehilohi are degp1 right hand side vectors of dimension dim;
 *   rhsrelolohi are degp1 right hand side vectors of dimension dim;
 *   rhsrehihilo are degp1 right hand side vectors of dimension dim;
 *   rhsrelohilo are degp1 right hand side vectors of dimension dim;
 *   rhsrehilolo are degp1 right hand side vectors of dimension dim;
 *   rhsrelololo are degp1 right hand side vectors of dimension dim;
 *   rhsimhihihi are degp1 right hand side vectors of dimension dim;
 *   rhsimlohihi are degp1 right hand side vectors of dimension dim;
 *   rhsimhilohi are degp1 right hand side vectors of dimension dim;
 *   rhsimlolohi are degp1 right hand side vectors of dimension dim;
 *   rhsimhihilo are degp1 right hand side vectors of dimension dim;
 *   rhsimlohilo are degp1 right hand side vectors of dimension dim;
 *   rhsimhilolo are degp1 right hand side vectors of dimension dim;
 *   rhsimlololo are degp1 right hand side vectors of dimension dim;
 *   solrehihihi are degp1 solution vectors of dimension dim;
 *   solrelohihi are degp1 solution vectors of dimension dim;
 *   solrehilohi are degp1 solution vectors of dimension dim;
 *   solrelolohi are degp1 solution vectors of dimension dim;
 *   solrehihilo are degp1 solution vectors of dimension dim;
 *   solrelohilo are degp1 solution vectors of dimension dim;
 *   solrehilolo are degp1 solution vectors of dimension dim;
 *   solrelololo are degp1 solution vectors of dimension dim;
 *   solimhihihi are degp1 solution vectors of dimension dim;
 *   solimlohihi are degp1 solution vectors of dimension dim;
 *   solimhilohi are degp1 solution vectors of dimension dim;
 *   solimlolohi are degp1 solution vectors of dimension dim;
 *   solimhihilo are degp1 solution vectors of dimension dim;
 *   solimlohilo are degp1 solution vectors of dimension dim;
 *   solimhilolo are degp1 solution vectors of dimension dim;
 *   solimlololo are degp1 solution vectors of dimension dim;
 *   resvecrehihihi has space for the residual power series;
 *   resvecrelohihi has space for the residual power series;
 *   resvecrehilohi has space for the residual power series;
 *   resvecrelolohi has space for the residual power series;
 *   resvecrehihilo has space for the residual power series;
 *   resvecrelohilo has space for the residual power series;
 *   resvecrehilolo has space for the residual power series;
 *   resvecrelololo has space for the residual power series;
 *   resvecimhihihi has space for the residual power series;
 *   resvecimlohihi has space for the residual power series;
 *   resvecimhilohi has space for the residual power series;
 *   resvecimlolohi has space for the residual power series;
 *   resvecimhihilo has space for the residual power series;
 *   resvecimlohilo has space for the residual power series;
 *   resvecimhilolo has space for the residual power series;
 *   resvecimlololo has space for the residual power series;
 *   totreslapsedms acculumates time spent by all kernels in milliseconds;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   resvecrehihihi are the highest doubles of the real parts
 *            of the residual series;
 *   resvecrelohihi are the second highest doubles of the real parts
 *            of the residual series;
 *   resvecrehilohi are the third highest doubles of the real parts
 *            of the residual series;
 *   resvecrelolohi are the fourth highest doubles of the real parts
 *            of the residual series;
 *   resvecrehihilo are the fourth lowest doubles of the real parts
 *            of the residual series;
 *   resvecrelohilo are the third lowest doubles of the real parts
 *            of the residual series;
 *   resvecrehilolo are the second lowest doubles of the real parts
 *            of the residual series;
 *   resvecrelololo are the lowest doubles of the real parts
 *            of the residual series;
 *   resvecimhihihi are the highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlohihi are the second highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimhilohi are the third highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlohihi are the fourth highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimhihilo are the fourth lowest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlohilo are the third lowest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimhilolo are the second lowest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlololo are the lowest doubles of the imaginary parts
 *            of the residual series;
 *   resmaxhihihi is the highest double of the max norm of the residual;
 *   resmaxlohihi is the second highest double of the max norm of the residual;
 *   resmaxhilohi is the third highest double of the max norm of the residual;
 *   resmaxlolohi is the fourth highest double of the max norm of the residual;
 *   resmaxhihilo is the fourth lowest double of the max norm of the residual;
 *   resmaxlohilo is the third lowest double of the max norm of the residual;
 *   resmaxhilolo is the second lowest double of the max norm of the residual;
 *   resmaxlololo is the lowest double of the max norm of the residual;
 *   totreslapsedms acculumates time spent by all kernels in milliseconds;
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

#endif
