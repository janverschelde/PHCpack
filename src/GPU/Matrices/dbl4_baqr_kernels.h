/* The file dbl4_baqr_kernels.h specifies functions for the
 * blocked accelerated QR in quad double precision. */

#ifndef __dbl4_baqr_kernels_h__
#define __dbl4_baqr_kernels_h__

/* The constants qd_shmemsize and cqd_shmemsize,
 * respectively for real and complex data, determine the upper bounds
 * on the size of the largest vectors the small kernels can handle.
 * The bounds were set experimentally, based on the available amount
 * of shared memory. */

#define qd_shmemsize 512
#define cqd_shmemsize 256

/* The constants inner_qd_shmemsize and outer_qd_shmemsize,
 * respectively for the many blocks of threads and the accumulating kernel,
 * determine how munch the shared memory is consumed. */

#define inner_qd_shmemsize 256
#define outer_qd_shmemsize 256

/* Both constants correspond to the number of threads in a block,
 * which typically are multiples of 32. */

__global__ void dbl4_small_house
 ( double *x0hihi, double *x0lohi, double *x0hilo, double *x0lolo,
   double *x1hihi, double *x1lohi, double *x1hilo, double *x1lolo,
   int dim, int dimLog2,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of a vector x of dimension dim+1,
 *   with one block of dim threads, on real data.
 *
 * ON ENTRY :
 *   x0hihi   highest double of the first element of the vector x; 
 *   x0lohi   second highest double of the first element of the vector x; 
 *   x0hilo   second lowest double of the first element of the vector x; 
 *   x0lolo   lowest double of the first element of the vector x; 
 *   x1hihi   array of dim doubles, with the highest doubles of x;
 *   x1lohi   array of dim doubles, with the second highest doubles of x;
 *   x1hilo   array of dim doubles, with the second lowest doubles of x;
 *   x1lolo   array of dim doubles, with the lowest doubles of x;
 *   dim      the dimension of the vector must equal the block size;
 *   dimLog2  equals ceil(log2((double) dim), used in sum reduction;
 *   vhihi    space allocated for dim+1 doubles;
 *   vlohi    space allocated for dim+1 doubles;
 *   vhilo    space allocated for dim+1 doubles;
 *   vlolo    space allocated for dim+1 doubles.
 *
 * ON RETURN :
 *   vhihi    highest doubles of the Householder vector;
 *   vlohi    second highest doubles of the Householder vector;
 *   vhilo    second lowest doubles of the Householder vector;
 *   vlolo    lowest doubles of the Householder vector;
 *   betahihi is the highest double of 2/(transpose(v)*v);
 *   betalohi is the second highest double of 2/(transpose(v)*v);
 *   betahilo is the second lowest double of 2/(transpose(v)*v);
 *   betalolo is the lowest double of 2/(transpose(v)*v). */

__global__ void cmplx4_small_house
 ( double *x0rehihi, double *x0relohi, double *x0rehilo, double *x0relolo,
   double *x0imhihi, double *x0imlohi, double *x0imhilo, double *x0imlolo,
   double *x1rehihi, double *x1relohi, double *x1rehilo, double *x1relolo,
   double *x1imhihi, double *x1imlohi, double *x1imhilo, double *x1imlolo,
   int dim, int dimLog2,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of a vector x of dimension dim+1,
 *   with one block of dim threads, on real data.
 *
 * ON ENTRY :
 *   x0rehihi highest double of the real part of the first element of x; 
 *   x0relohi second highest double of the real part of the first element of x; 
 *   x0rehilo second lowest double of the real part of the first element of x; 
 *   x0relolo lowest double of the real part of the first element of x; 
 *   x0imhihi highest double of the imaginary parts of the first element of x; 
 *   x0imlohi second highest double of the imaginary parts
 *            of the first element of x; 
 *   x0imhilo second lowest double of the imaginary parts
 *            of the first element of x; 
 *   x0imlolo lowest double of the imaginary parts of the first element of x; 
 *   x1rehihi array of dim highest doubles of the real parts of x;
 *   x1relohi array of dim second highest doubles of the real parts of x;
 *   x1rehilo array of dim second lowest doubles of the real parts of x;
 *   x1relolo array of dim lowest doubles of the real parts of x;
 *   x1imhihi array of dim highest doubles of the imaginary parts of x;
 *   x1imlohi array of dim second highest doubles of the imaginary parts of x;
 *   x1imhilo array of dim second lowest doubles of the imaginary parts of x;
 *   x1imlolo array of dim lowest doubles of the imaginary parts of x;
 *   dim      the dimension of the vector must equal the block size;
 *   dimLog2  equals ceil(log2((double) dim), used in sum reduction;
 *   vrehihi  space allocated for dim+1 doubles;
 *   vrelohi  space allocated for dim+1 doubles;
 *   vrehilo  space allocated for dim+1 doubles;
 *   vrelolo  space allocated for dim+1 doubles;
 *   vimhihi  space allocated for dim+1 doubles;
 *   vimlohi  space allocated for dim+1 doubles;
 *   vimhilo  space allocated for dim+1 doubles;
 *   vimlolo  space allocated for dim+1 doubles.
 *
 * ON RETURN :
 *   vrehihi  highest doubles of the real parts of the Householder vector;
 *   vrelohi  second highest doubles of the real parts
 *            of the Householder vector;
 *   vrehilo  second lowest doubles of the real parts
 *            of the Householder vector;
 *   vrelolo  lowest doubles of the real parts of the Householder vector;
 *   vimhihi  highest doubles of the imaginary parts of the Householder vector;
 *   vimlohi  second highest doubles of the imaginary parts
 *            of the Householder vector;
 *   vimhilo  second lowest doubles of the imaginary parts
 *            of the Householder vector;
 *   vimlolo  lowest doubles of the imaginary parts of the Householder vector;
 *   betahihi is the highest double of 2/(transpose(v)*v);
 *   betalohi is the second highest double of 2/(transpose(v)*v);
 *   betahilo is the second lowest double of 2/(transpose(v)*v);
 *   betalolo is the lowest double of 2/(transpose(v)*v). */

__global__ void dbl4_large_sum_of_squares
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int dim, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in a vector,
 *   given with high doubles in vhi and low doubles in vlo,
 *   as needed in the 2-norm of the vector, with many blocks,
 *   on real data.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vhihi     highest doubles of a vector v;
 *   vlohi     second highest doubles of a vector v;
 *   vhilo     second lowest doubles of v;
 *   vlolo     lowest doubles of v;
 *   sumshihi  space for as many doubles as the number of blocks;
 *   sumslohi  space for as many doubles as the number of blocks;
 *   sumshilo  space for as many doubles as the number of blocks;
 *   sumslolo  space for as many doubles as the number of blocks;
 *   dim       number of elements in v;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshihi  highest doubles of computed sums of squares of vector slices,
 *             the i-th entry is computed by the i-th block of threads;
 *   sumslohi  second highest doubles of computed sums of squares;
 *   sumshilo  second lowest doubles of computed sums of squares;
 *   sumslolo  lowest doubles of computed sums of squares. */

__global__ void cmplx4_large_sum_of_squares
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int dim, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in a vector,
 *   given with high doubles in vhi and low doubles in vlo,
 *   as needed in the 2-norm of the vector, with many blocks,
 *   on complex data.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vrehihi   highest doubles of the real parts of a vector v;
 *   vrelohi   second highest doubles of the real parts of a vector v;
 *   vrehilo   second lowest doubles of the real parts of v;
 *   vrelolo   lowest doubles of the real parts of v;
 *   vimhihi   highest doubles of the imaginary parts of v;
 *   vimlohi   second highest doubles of the imaginary parts of v;
 *   vimhilo   second lowest doubles of the imaginary parts of v;
 *   vimlolo   lowest doubles of the imaginary parts of v;
 *   sumshihi  has space for as many doubles as the number of blocks;
 *   sumslohi  has space for as many doubles as the number of blocks;
 *   sumslohi  has space for as many doubles as the number of blocks;
 *   sumslolo  has space for as many doubles as the number of blocks;
 *   dim       number of elements in v;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshihi  highest doubles of computed sums of squares of vector slices,
 *             the i-th entry is computed by the i-th block of threads;
 *   sumslohi  second highest doubles of computed sums of squares;
 *   sumshilo  second lowest doubles of computed sums of squares;
 *   sumslolo  lowest doubles of computed sums of squares. */

__global__ void dbl4_sum_accumulator
 ( double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int nbsums, int nbsumsLog2,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo );
/*
 * DESCRIPTION :
 *   Accumulates the sum by one block of threads, on real data.
 *
 * REQUIRED : the number of threads should be equal to nbsums.
 *
 * ON ENTRY :
 *   sumshihi  highest doubles of computed sums of squares of vector slices;
 *   sumslohi  second highest doubles of computed sums of squares;
 *   sumshilo  second lowest doubles of computed sums of squares;
 *   sumslolo  lowest doubles of computed sums of squares;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction.
 *
 * ON RETURN :
 *   acchihi   the highest double of the sum;
 *   acclohi   the second highest double of the sum;
 *   acchilo   the second lowest double of the sum;
 *   acclolo   the lowest double of the sum. */

__global__ void dbl4_normalize
 ( int dim, int szt,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *v0hihi, double *v0lohi, double *v0hilo, double *v0lolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo );
/*
 * DESCRIPTION :
 *   Divides every element in the vector x with the same number v0,
 *   using multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   dim       number of elements in the vectors x and v;
 *   szt       size of one block;
 *   xhihi     highest doubles of the vector x;
 *   xlohi     second highest doubles of the vector x;
 *   xhilo     second lowest doubles of the vector x;
 *   xlolo     lowest doubles of the vector x;
 *   v0hihi    highest double of v0;
 *   v0lohi    second highest double of v0;
 *   v0hilo    second lowest double of v0;
 *   v0lolo    lowest double of v0;
 *   vhihi     space for dim doubles;
 *   vlohi     space for dim doubles;
 *   vhilo     space of dim doubles;
 *   vlolo     space of dim doubles.
 *
 * ON RETURN :
 *   vhihi     highest doubles of x divided by v0;
 *   vlohi     second highest doubles of x divided by v0;
 *   vhilo     second lowest doubles of x divided by v0;
 *   vlolo     lowest doubles of x divided by v0. */

__global__ void cmplx4_normalize
 ( int dim, int szt,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *inv0rehihi, double *inv0relohi,
   double *inv0rehilo, double *inv0relolo,
   double *inv0imhihi, double *inv0imlohi,
   double *inv0imhilo, double *inv0imlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo );
/*
 * DESCRIPTION :
 *   Multiplies every element in the vector x with the same number v0,
 *   using multiple blocks of threads, on complex data.
 *
 * ON ENTRY :
 *   dim       number of elements in the vectors x and v;
 *   szt       size of one block;
 *   xrehihi   highest doubles of the real parts of x;
 *   xrelohi   second highest doubles of the real parts of x;
 *   xrehilo   second lowest doubles of the real parts of x;
 *   xrelolo   lowest doubles of the real parts of x;
 *   ximhihi   highest doubles of the imaginary parts of x;
 *   ximlohi   second highest doubles of the imaginary parts of x;
 *   ximhilo   second lowest doubles of the imaginary parts of x;
 *   ximlolo   lowest doubles of the imaginary parts of x;
 *   inv0rehihi is the highest double of the real part of 1/v0;
 *   inv0relohi is the second highest double of the real part of 1/v0;
 *   inv0rehilo is the second lowest double of the real part of 1/v0;
 *   inv0relolo is the lowest double of the real part of 1/v0;
 *   inv0imhihi is the highest double of the imaginary part of 1/v0;
 *   inv0imlohi is the second highest double of the imaginary part of 1/v0;
 *   inv0imhilo is the second lowest double of the imaginary part of 1/v0;
 *   inv0imlolo is the lowest double of the imaginary part of 1/v0;
 *   vrehihi   space for dim doubles;
 *   vrelohi   space for dim doubles;
 *   vrehilo   space for dim doubles;
 *   vrelolo   space for dim doubles;
 *   vimhihi   space for dim doubles;
 *   vimlohi   space for dim doubles;
 *   vimhilo   space for dim doubles;
 *   vimlolo   space for dim doubles.
 *
 * ON RETURN :
 *   vrehihi   highest doubles of x multiplied by 1/v0;
 *   vrelohi   second highest doubles of x multiplied by 1/v0;
 *   vrehilo   second lowest doubles of x multiplied by 1/v0;
 *   vrelolo   lowest doubles of x multiplied by 1/v0;
 *   vimhihi   highest doubles of x multiplied by 1/v0;
 *   vimlohi   second highest doubles of x multiplied by 1/v0;
 *   vimhilo   second lowest doubles of x multiplied by 1/v0;
 *   vimlolo   lowest doubles of x multiplied by 1/v0. */

__global__ void dbl4_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi,
   double *betahilo, double *betalolo );
/*
 * DESCRIPTION :
 *   Updates the matrix R starting at column k
 *   with the Householder vector in v and beta, on real data,
 *   with a total of ncols - k threads in one block.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   Rhihi    highest doubles of R, stored column wise;
 *   Rlohi    second highest doubles of R, stored column wise;
 *   Rhilo    second lowest doubles of R, stored column wise;
 *   Rlolo    lowest doubles of R, stored column wise;
 *   vhihi    highest doubles of the Householder vector;
 *   vlohi    second highest doubles of the Householder vector;
 *   vhilo    second lowest doubles of the Householder vector;
 *   vlolo    lowest doubles of the Householder vector;
 *   betahihi is the highest double of 2/(transpose(v)*v);
 *   betalohi is the second highest double of 2/(transpose(v)*v);
 *   betahilo is the second lowest double of 2/(transpose(v)*v).
 *   betalolo is the lowest double of 2/(transpose(v)*v).
 *
 * ON RETURN :
 *   Rhihi    highest doubles of the updated R, which is trapezoidal;
 *   Rlohi    second highest doubles of the updated R;
 *   Rhilo    second lowest doubles of the updated R;
 *   Rlolo    lowest doubles of the updated R. */

__global__ void cmplx4_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo );
/*
 * DESCRIPTION :
 *   Updates the matrix R starting at column k
 *   with the Householder vector in v and beta, on complex data,
 *   with a total of ncols - k threads in one block.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   Rrehihi  highest doubles of the real parts of R, stored column wise;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R.
 *
 * ON RETURN :
 *   Rrehihi  highest doubles of the real parts of the updated R;
 *   Rrelohi  second highest doubles of the real parts of the updated R;
 *   Rrehilo  second lowest doubles of the real parts of the updated R;
 *   Rrelolo  lowest doubles of the real parts of the updated R;
 *   Rimhihi  highest doubles of the imaginary parts of the updated R;
 *   Rimlohi  second highest doubles of the imaginary parts of the updated R;
 *   Rimhilo  second lowest doubles of the imaginary parts of the updated R;
 *   Rimlolo  lowest doubles of the imaginary parts of the updated R. */

__global__ void dbl4_RTdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *RTdotvhihi, double *RTdotvlohi,
   double *RTdotvhilo, double *RTdotvlolo );
/*
 * DESCRIPTION :
 *   The elements of the matrix RTdotv are the elements of R^T,
 *   multiplied with the corresponding element of v.
 *   Multiple blocks of threads operate on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   szt      size of one tile and number of threads in one block;
 *   colidx   index of the current column in R;
 *   Roffset  offset in R for the first row to start;
 *   dim      number of columns in R;
 *   Rhihi    highest doubles of acolumn wise stored matrix 
 *            with number of rows equal to nrows;
 *   Rlohi    second highest doubles of acolumn wise stored matrix 
 *            with number of rows equal to nrows;
 *   Rhilo    second lowest doubles of acolumn wise stored matrix 
 *            with number of rows equal to nrows;
 *   Rlolo    lowest doubles of acolumn wise stored matrix 
 *            with number of rows equal to nrows;
 *   vhihi    start of the highest doubles of the first nonzero element
 *            of a Householder vector;
 *   vlohi    start of the second highest doubles of the first nonzero element
 *            of a Householder vector;
 *   vhilo    start of the second lowest doubles of the first nonzero element
 *            of a Householder vector;
 *   vlolo    start of the lowest doubles of the first nonzero element
 *            of a Householder vector;
 *   RTdotvhihi has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvlohi has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvhilo has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvlolo has space for a matrix of nrows-by-szt, plus some padding.
 *
 * ON RETURN :
 *   RTdotvhihi are the highest doubles of the element-by-element products
 *            of R^T with v, stored row by row;
 *   RTdotvlohi are the second highest doubles of the element-by-element
 *            products of R^T with v, stored row by row;
 *   RTdotvhilo are the second lowest doubles of the element-by-element
 *            products of R^T with v, stored row by row;
 *   RTdotvlolo are the lowest doubles of the element-by-element products
 *            of R^T with v, stored row by row. */

__global__ void cmplx4_RHdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *RHdotvrehihi, double *RHdotvrelohi,
   double *RHdotvrehilo, double *RHdotvrelolo,
   double *RHdotvimhihi, double *RHdotvimlohi,
   double *RHdotvimhilo, double *RHdotvimlolo );
/*
 * DESCRIPTION :
 *   The elements of the matrix RHdotv are the elements of R^H,
 *   multiplied with the corresponding element of v.
 *   Multiple blocks of threads operate on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   szt      size of one tile and number of threads in one block;
 *   colidx   index of the current column in R;
 *   Roffset  offset in R for the first row to start;
 *   dim      number of columns in R;
 *   Rrehihi  highest doubles of the real parts of R, a column wise stored
 *            matrix with number of rows equal to nrows;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R;
 *   vrehihi  start of the the highest doubles of the real parts of the
 *            first nonzero element of a Householder vector v;
 *   vrelohi  second highest doubles of the real parts of v;
 *   vrehilo  second lowest doubles of the real parts of v;
 *   vrelolo  lowest doubles of the real parts of v;
 *   vimhihi  highest doubles of the imaginary parts of v;
 *   vimlohi  second highest doubles of the imaginary parts of v;
 *   vimhilo  second lowest doubles of the imaginary parts of v;
 *   vimlolo  lowest doubles of the imaginary parts of v;
 *   RHdotvrehihi has space for a matrix of nrows-by-szt, plus some padding;
 *   RHdotvrelohi has space for a matrix of nrows-by-szt, plus some padding;
 *   RHdotvrehilo has space for a matrix of nrows-by-szt, plus some padding;
 *   RHdotvrelolo has space for a matrix of nrows-by-szt, plus some padding;
 *   RHdotvimhihi has space for a matrix of nrows-by-szt, plus some padding;
 *   RHdotvimlohi has space for a matrix of nrows-by-szt, plus some padding;
 *   RHdotvimhilo has space for a matrix of nrows-by-szt, plus some padding.
 *   RHdotvimlolo has space for a matrix of nrows-by-szt, plus some padding.
 *
 * ON RETURN :
 *   RHdotvrehihi has the highest doubles of the real parts of the 
 *            element-by-element products of R^H with v, stored row by row;
 *   RHdotvrelohi has the second highest doubles of the real parts of RHdotv;
 *   RHdotvrehilo has the second lowest doubles of the real parts of RHdotv;
 *   RHdotvrelolo has the lowest doubles of the real parts of RHdotv;
 *   RHdotvimhihi has the highest doubles of the imag parts of RHdotv;
 *   RHdotvimlohi has the second highest doubles of the imag parts of RHdotv;
 *   RHdotvimhilo has the second lowest doubles of the imag parts of RHdotv;
 *   RHdotvimlolo has the lowest doubles of the imag parts of RHdotv. */

__global__ void cmplx4_RHdotvre
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *RHdotvrehihi, double *RHdotvrelohi,
   double *RHdotvrehilo, double *RHdotvrelolo,
   double *RHdotvimhihi, double *RHdotvimlohi,
   double *RHdotvimhilo, double *RHdotvimlolo );
/*
 * DESCRIPTION :
 *   Computes only the real parts of R^H dotted with v. */

__global__ void cmplx4_RHdotvim
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *RHdotvimhihi, double *RHdotvimlohi,
   double *RHdotvimhilo, double *RHdotvimlolo );
/*
 * DESCRIPTION :
 *   Computes only the imaginary parts of R^H dotted with v. */

__global__ void dbl4_sum_betaRTdotv
 ( int nrows,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *RTdotvhihi, double *RTdotvlohi,
   double *RTdotvhilo, double *RTdotvlolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo );
/*
 * DESCRIPTION :
 *   Adds the rows in RTdotv to obtain w = beta*R^T*v,
 *   with one block of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in RTdotv;
 *   betahihi is the highest double of the beta for the Householder vector;
 *   betalohi is the second highest double of the beta;
 *   betahilo is the second lowest double of the beta;
 *   betahilo is the lowest double of the beta;
 *   RTdotvhihi has the highest doubles of all products of the elements
 *            of R dotted v;
 *   RTdotvlohi has the second highest doubles of of R dotted v;
 *   RTdotvhilo has the scond lowest doubles of of R dotted v;
 *   RTdotvlolo has the lowest doubles of of R dotted v;
 *   whihi    space for the highest doubles of beta*R^T*v.
 *   wlohi    space for the second highest doubles of beta*R^T*v.
 *   whilo    space for the second lowest doubles of beta*R^T*v.
 *   wlolo    space for the lowest doubles of beta*R^T*v.
 *
 * ON RETURN :
 *   whihi    the highest doubles of beta*R^T*v;
 *   wlohi    the second highest doubles of beta*R^T*v;
 *   whilo    the second lowest doubles of beta*R^T*v;
 *   wlolo    the lowest doubles of beta*R^T*v. */

__global__ void cmplx4_sum_betaRHdotv
 ( int nrows,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *RTdotvrehihi, double *RTdotvrelohi,
   double *RTdotvrehilo, double *RTdotvrelolo,
   double *RTdotvimhihi, double *RTdotvimlohi,
   double *RTdotvimhilo, double *RTdotvimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo );
/*
 * DESCRIPTION :
 *   Adds the rows in RHdotv to obtain w = beta*R^H*v,
 *   with one block of threads, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in RHdotv;
 *   betahihi is the highest double of the beta for the Householder vector;
 *   betalohi is the second highest double of the beta;
 *   betahilo is the second lowest double of the beta;
 *   betahilo is the lowest double of the beta;
 *   RTdotvrehihi has the highest doubles of the real parts of all products
 *            of the elements of R dotted v;
 *   RTdotvrelohi has the second highest doubles of the real parts
 *            of R dotted v;
 *   RTdotvrehilo has the second lowest doubles of the real parts
 *            of R dotted v;
 *   RTdotvrelolo has the lowest doubles of the real parts of R dotted v;
 *   RTdotvimhihi has the highest doubles of the imaginary parts
 *            of R dotted v;
 *   RTdotvimhihi has the second highest doubles of the imaginary parts
 *            of R dotted v;
 *   RTdotvimhilo has the second lowest doubles of the imaginary parts
 *            of R dotted v;
 *   RTdotvimlolo has the lowest doubles of the imaginary parts
 *            of R dotted v;
 *   wrehihi  space for the highest doubles of the real parts of beta*R^H*v;
 *   wrelohi  space for the second highest doubles of the real parts
 *            of beta*R^H*v;
 *   wrehilo  space for the second lowest doubles of the real parts
 *            of beta*R^H*v;
 *   wrelolo  space for the lowest doubles of the real parts of beta*R^H*v;
 *   wimhihi  space for the highest doubles of the imaginary parts
 *            of beta*R^H*v;
 *   wimlohi  space for the second highest doubles of the imaginary parts
 *            of beta*R^H*v;
 *   wimhilo  space for the second lowest doubles of the imaginary parts
 *            of beta*R^H*v.
 *   wimlolo  space for the lowest doubles of the imaginary parts
 *            of beta*R^H*v.
 *
 * ON RETURN :
 *   wrehihi  highest doubles of the real parts of beta*R^H*v;
 *   wrelohi  second highest doubles of the real parts of beta*R^H*v;
 *   wrehilo  second lowest doubles of the real parts of beta*R^H*v;
 *   wrelolo  lowest doubles of the real parts of beta*R^H*v;
 *   wimhihi  highest doubles of the imaginary parts of beta*R^H*v;
 *   wimlohi  second highest doubles of the imaginary parts of beta*R^H*v;
 *   wimhilo  second lowest doubles of the imaginary parts of beta*R^H*v;
 *   wimlolo  lowest doubles of the imaginary parts of beta*R^H*v. */

__global__ void dbl4_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo );
/*
 * DESCRIPTION :
 *   Applies the Householder vector to update R, on real data.
 *
 * REQUIRED : nrows - k > szt as multiple blocks are used.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   Rhihi    highest doubles of the nrows-by-ncols matrix R,
 *            stored column wise;
 *   Rlohi    second highest doubles of R;
 *   Rhilo    second lowest doubles of R;
 *   Rlolo    lowest doubles of R;
 *   vhihi    highest double of the Householder vector;
 *   vlohi    second highest double of the Householder vector;
 *   vhilo    second lowest double of the Householder vector;
 *   vlolo    lowest double of the Householder vector;
 *   betahihi has the highest double of the beta corresponding with v;
 *   betalohi has the second highest double of the beta;
 *   betalohi has the second lowest double of the beta;
 *   betalolo has the lowest double of the beta;
 *   whihi    highest doubles of beta*R^T*v;
 *   wlohi    second highest doubles of beta*R^T*v;
 *   whilo    second lowest doubles of beta*R^T*v.
 *   wlolo    lowest doubles of beta*R^T*v.
 *
 * ON RETURN :
 *   Rhihi    highest doubles of the updated R;
 *   Rlohi    second highest doubles of the updated R;
 *   Rhilo    second lowest doubles of the updated R;
 *   Rlolo    lowest doubles of the updated R. */

__global__ void cmplx4_medium_subvbetaRHv
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo );
/*
 * DESCRIPTION :
 *   Applies the Householder vector to update R, on complex data.
 *
 * REQUIRED : nrows - k > szt as multiple blocks are used.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   Rrehihi  highest doubles of the real parts of R,
 *            an nrows-by-ncols matrix, stored column wise;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R;
 *   vrehi    high doubles of the real parts of the Householder vector;
 *   vrelo    low doubles of the real parts of the Householder vector;
 *   vimhi    high doubles of the imaginary parts of the Householder vector;
 *   vimlo    low doubles of the imaginary parts of the Householder vector;
 *   betahi   high double of the beta corresponding with v;
 *   betalo   low double of the beta corresponding with v;
 *   wrehi    high doubles of the real parts of beta*R^H*v;
 *   wrelo    low doubles of the real parts of beta*R^H*v;
 *   wimhi    high doubles of the imaginary parts of beta*R^H*v;
 *   wimlo    low doubles of the imaginary parts of beta*R^H*v.
 *
 * ON RETURN :
 *   Rrehihi  highest doubles of the real parts of the updated R;
 *   Rrelohi  second highest doubles of the real parts of the updated R;
 *   Rrehilo  second lowest doubles of the real parts of the updated R;
 *   Rrelolo  lowest doubles of the real parts of the updated R;
 *   Rimhihi  highest doubles of the imaginary parts of the updated R;
 *   Rimlohi  second highest doubles of the imaginary parts of the updated R;
 *   Rimhilo  second lowest doubles of the imaginary parts of the updated R;
 *   Rimlolo  lowest doubles of the imaginary parts of the updated R. */

__global__ void cmplx4_medium_subvbetaRHvRe
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo );
/*
 * DESCRIPTION :
 *   Updates the real parts of R. */

__global__ void cmplx4_medium_subvbetaRHvIm
 ( int nrows, int ncols, int szt, int k,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo );
/*
 * DESCRIPTION :
 *    Updates the imaginary parts of R. */

__global__ void dbl4_beta_times_V
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo );
/*
 * DESCRIPTION :
 *   Computes the first vector in the W representation of the Householder
 *   transformations, multiplying B[0] with the first vector of V,
 *   and flipping the sign, with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V and W;
 *   szt      size of one tile and the number of threads in a block;
 *   Bhihi    highest double of the betas computed by house;
 *   Blohi    second highest double of the betas;
 *   Bhilo    second lowest double of the betas;
 *   Blolo    lowest double of the betas;
 *   Vhihi    Vhihi[nrows*i] is the start of the highest doubles of the i-th
 *            Householder vector, with i zeros inserted so V is trapezoidal;
 *   Vlohi    second highest doubles of V;
 *   Vhilo    second lowest doubles of V;
 *   Vlolo    lowest doubles of V;
 *   Whihi    space for the highest doubles of the W matrix;
 *   Wlohi    space for the highest doubles of the W matrix;
 *   Whilo    space for the lowest doubles of the W matrix;
 *   Wlolo    space for the lowest doubles of the W matrix.
 *
 * ON RETURN :
 *   Whi      the first nrows numbers store the highest doubles of the
 *            first vector of the W matrix in the WY representation;
 *   Whilo    the second highest doubles of W;
 *   Whilo    the second lowest doubles of W;
 *   Wlolo    the lowest doubles of W. */

__global__ void cmplx4_beta_times_V
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo );
/*
 * DESCRIPTION :
 *   Computes the first vector in the W representation of the Householder
 *   transformations, multiplying B[0] with the first vector of V,
 *   and flipping the sign, with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V and W;
 *   szt      size of one tile and the number of threads in a block;
 *   Bhihi    highest double of the betas computed by house;
 *   Blohi    second highest double of the betas;
 *   Bhilo    second lowest double of the betas;
 *   Blolo    lowest double of the betas;
 *   Vrehihi  Vrehihi[nrows*i] is the start of the high doubles of the real
 *            parts of the i-th Householder vector, with i zeros inserted;
 *   Vrelohi  second highest doubles of the real parts of V;
 *   Vrehilo  second lowest doubles of the real parts of V;
 *   Vrelolo  lowest doubles of the real parts of V;
 *   Vimhihi  highest doubles of the imaginary parts of V;
 *   Vimlohi  second highest doubles of the imaginary parts of V;
 *   Vimhilo  second lowest doubles of the imaginary parts of V;
 *   Vimlolo  lowest doubles of the imaginary parts of V;
 *   Wrehihi  space for the highest doubles of the real parts of the W;
 *   Wrelohi  space for the second highest doubles
 *            of the real parts of the W;
 *   Wrehilo  space for the second lowest doubles
 *            of the real parts of the W;
 *   Wrelolo  space for the lowest doubles of the real parts of the W;
 *   Wimhihi  space for the highest doubles of the imaginary parts of W;
 *   Wimlohi  space for the second highest doubles
 *            of the imaginary parts of W;
 *   Wimhilo  space for the second lowest doubles
 *            of the imaginary parts of W.
 *   Wimlolo  space for the lowest doubles of the imaginary parts of W.
 *
 * ON RETURN :
 *   Wrehihi  the first nrows numbers store the highest doubles of the real
 *            parts of the first vector of the W the WY representation;
 *   Wrelohi  second highest doubles of the real parts of W;
 *   Wrehilo  second lowest doubles of the real parts of W;
 *   Wrelolo  lowest doubles of the real parts of W;
 *   Wimhihi  highest doubles of the imaginary parts of W;
 *   Wimlohi  second highest doubles of the imaginary parts of W;
 *   Wimhilo  second lowest doubles of the imaginary parts of W;
 *   Wimlolo  lowest doubles of the imaginary parts of W. */

__global__ void dbl4_initialize_WYT
 ( int dim, int szt,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo );
/*
 * DESCRIPTION :
 *   Initializes the matrix YWT with the product of the first dim products
 *   of W with Y elements with multiple blocks of szt threads,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in WYT;
 *   szt      number of threads in one block;
 *   Vhihi    the first dim numbers define the highest doubles of V,
 *            the first Householder vector;
 *   Vlohi    the first dim numbers define the second highest doubles of V;
 *   Vhilo    the first dim numbers define the second lowest doubles of V;
 *   Vlolo    the first dim numbers define the lowest doubles of V;
 *   Whihi    the first dim numbers define the highest doubles of
 *            the first column in the W matrix;
 *   Wlohi    the first dim numbers define the second highest doubles of W;
 *   Whilo    the first dim numbers define the second lowest doubles of W;
 *   Wlolo    the first dim numbers define the lowest doubles of W;
 *   WYThihi  space for a dim-by-dim matrix;
 *   WYTlohi  space for a dim-by-dim matrix;
 *   WYThilo  space for a dim-by-dim matrix;
 *   WYTlolo  space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   WYThihi  highest doubles of w*y^T,
 *            where y and w are the first columns of V and W;
 *   WYTlohi  second highest doubles of w*y^T;
 *   WYThilo  second lowest doubles of w*y^T;
 *   WYTlolo  lowest doubles of w*y^T. */

__global__ void cmplx4_initialize_WYH
 ( int dim, int szt,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYTrehihi, double *WYTrelohi,
   double *WYTrehilo, double *WYTrelolo,
   double *WYTimhihi, double *WYTimlohi,
   double *WYTimhilo, double *WYTimlolo );
/*
 * DESCRIPTION :
 *   Initializes the matrix YWT with the product of the first dim products
 *   of W with Y elements with multiple blocks of szt threads,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in WYT;
 *   szt      number of threads in one block;
 *   Vrehihi  the first dim numbers define the highest doubles of
 *            the real parts of V, the first Householder vector;
 *   Vrelohi  second highest doubles of the real parts of V;
 *   Vrehilo  second lowest doubles of the real parts of V;
 *   Vrelolo  lowest doubles of the real parts of V;
 *   Vimhihi  highest doubles of the imaginary parts of V;
 *   Vimlohi  second highest doubles of the imaginary parts of V;
 *   Vimhilo  second lowest doubles of the imaginary parts of V;
 *   Vimlolo  lowest doubles of the imaginary parts of V;
 *   Wrehihi  the first dim numbers define the highest doubles of
 *            the real parts of the first column in the W matrix;
 *   Wrelohi  second highest doubles of the real parts of W;
 *   Wrehilo  second lowest doubles of the real parts of W;
 *   Wrelolo  lowest doubles of the real parts of W;
 *   Wimhihi  highest doubles of the imaginary parts of W;
 *   Wimlohi  second highest doubles of the imaginary parts of W;
 *   Wimhilo  second lowest doubles of the imaginary parts of W;
 *   Wimlolo  lowest doubles of the imaginary parts of W;
 *   WYTrehihi has space for a dim-by-dim matrix;
 *   WYTrelohi has space for a dim-by-dim matrix;
 *   WYTrehilo has space for a dim-by-dim matrix;
 *   WYTrelolo has space for a dim-by-dim matrix;
 *   WYTimhihi has space for a dim-by-dim matrix;
 *   WYTimlohi has space for a dim-by-dim matrix;
 *   WYTimhilo has space for a dim-by-dim matrix;
 *   WYTimlolo has space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   WYTrehihi has the highest doubles of the real part of w*y^T,
 *            where y and w are the first columns of V and W;
 *   WYTrelohi has the second highest doubles of the real part of w*y^T;
 *   WYTrehilo has the second lowest doubles of the real part of w*y^T;
 *   WYTrelolo has the lowest doubles of the real part of w*y^T;
 *   WYTimhihi has the highest doubles of the imaginary part of w*y^T;
 *   WYTimlohi has the second highest doubles of the imaginary part of w*y^T;
 *   WYTimhilo has the second lowest doubles of the imaginary part of w*y^T;
 *   WYTimlolo has the lowest doubles of the imaginary part of w*y^T. */

__global__ void dbl4_update_WYT
 ( int dim, int szt,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo );
/*
 * DESCRIPTION :
 *   Updates the matrix WYT with the product of the first dim products
 *   of W with Y elements with multiple blocks of szt threads,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in YWT;
 *   szt      number of threads in one block;
 *   Vhihi    the first dim numbers define the highest doubles
 *            of the next Householder vector;
 *   Vlohi    second highest doubles of V;
 *   Vhilo    second lowest doubles of V;
 *   Vlolo    lowest doubles of V;
 *   Whihi    the first dim numbers define the highest doubles
 *            of the next column in the W matrix;
 *   Wlohi    second highest doubles of W;
 *   Whilo    second lowest doubles of W;
 *   Wlolo    lowest doubles of W;
 *   WYThihi  dim-by-dim matrix with the highest doubles for WYT;
 *   WYTlohi  dim-by-dim matrix with the second highest doubles for WYT;
 *   WYThilo  dim-by-dim matrix with the second lowest doubles for WYT;
 *   WYTlolo  dim-by-dim matrix with the lowest doubles for WYT.
 *
 * ON RETURN :
 *   WYThihi  highest doubles of the updated matrix W*Y^T;
 *   WYThihi  second highest doubles of the updated matrix W*Y^T;
 *   WYTlolo  second lowest doubles of the updated matrix W*Y^T;
 *   WYTlolo  lowest doubles of the updated matrix W*Y^T. */

__global__ void cmplx4_update_WYH
 ( int dim, int szt,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYHrehihi, double *WYHrelohi,
   double *WYHrehilo, double *WYHrelolo,
   double *WYHimhihi, double *WYHimlohi,
   double *WYHimhilo, double *WYHimlolo );
/*
 * DESCRIPTION :
 *   Updates the matrix YWT with the product of the first dim products
 *   of W with Y elements with multiple blocks of szt threads,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in WYT;
 *   szt      number of threads in one block;
 *   Vrehi    the first dim numbers define the high doubles of
 *            the real parts of the first Householder vector;
 *   Vrelo    the first dim numbers define the low doubles of
 *            the real parts of the first Householder vector;
 *   Vimhi    the first dim numbers define the high doubles of
 *            the imaginary parts of the first Householder vector;
 *   Vimlo    the first dim numbers define the low doubles of
 *            the imaginary parts of the first Householder vector;
 *   Wrehi    the first dim numbers define the high doubles of
 *            the real parts of the first column in the W matrix;
 *   Wrelo    the first dim numbers define the low doubles of
 *            the real parts of the first column in the W matrix;
 *   Wimhi    the first dim numbers define the high doubles of
 *            the imaginary parts of the first column in the W matrix;
 *   Wimlo    the first dim numbers define the low doubles of
 *            the imaginary parts of the first column in the W matrix;
 *   WYHrehi  dim-by-dim matrix with the high doubles of the real
 *            parts of W*Y^H;
 *   WYHrelo  dim-by-dim matrix with the low doubles of the imaginary
 *            parts of W*Y^H;
 *   WYHimhi  dim-by-dim matrix with the high doubles of the real
 *            parts of W*Y^H;
 *   WYHimlo  dim-by-dim matrix with the low doubles of the imaginary
 *            parts of W*Y^H.
 *
 * ON RETURN :
 *   WYHrehi  high doubles of the real parts of the updated W*Y^H,
 *            where y and w are the first columns of V and W;
 *   WYHrelo  low doubles of the real parts of the updated W*Y^H,
 *            where y and w are the first columns of V and W;
 *   WYHimhi  high doubles of the imaginary parts of the updated W*Y^H,
 *            where y and w are the first columns of V and W;
 *   WYHimlo  low doubles of the imaginary part of the updated W*Y^H,
 *            where y and w are the first columns of V and W. */

__global__ void dbl4_beta_next_W
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo );
/*
 * DECRIPTION :
 *   Computes the next column in the W matrix, with multiple blocks,
 *   on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V, W, and the dimension of WYT;
 *   szt      number of threads in one block;
 *   Bhihi    highest double of the beta with the next Householder vector;
 *   Blohi    second highest double of the beta;
 *   Bhilo    second lowest double of the beta;
 *   Blolo    lowest double of the beta;
 *   Vhihi    the first dim numbers define the highest doubles of V,
 *            the next Householder vector;
 *   Vlohi    the first dim numbers define the second highest doubles of V;
 *   Vhilo    the first dim numbers define the second lowest doubles of V;
 *   Vlolo    the first dim numbers define the lowest doubles of V;
 *   Whihi    the first dim numbers define the highest doubles of
 *            the next column in the W matrix;
 *   Wlohi    the first dim numbers define the second highest doubles of W;
 *   Whilo    the first dim numbers define the second lowest doubles of W;
 *   Wlolo    the first dim numbers define the lowest doubles of W;
 *   WYThihi  dim-by-dim matrix with the highest doubles for WYT;
 *   WYTlohi  dim-by-dim matrix with the second highest doubles for WYT;
 *   WYThilo  dim-by-dim matrix with the second lowest doubles for WYT.
 *   WYTlolo  dim-by-dim matrix with the lowest doubles for WYT.
 *
 * ON RETURN :
 *   Whihi    the highest doubles of the next column of W;
 *   Wlohi    the second highest doubles of the next column of W;
 *   Whilo    the second lowest doubles of the next column of W;
 *   Wlolo    the lowest doubles of the next column of W. */

__global__ void cmplx4_beta_next_W
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYHrehihi, double *WYHrelohi,
   double *WYHrehilo, double *WYHrelolo,
   double *WYHimhihi, double *WYHimlohi,
   double *WYHimhilo, double *WYHimlolo );
/*
 * DECRIPTION :
 *   Computes the next column in the W matrix, with multiple blocks,
 *   on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V, W, and the dimension of WYT;
 *   szt      number of threads in one block;
 *   Bhihi    highest double of the beta with the next Householder vector;
 *   Blohi    second highest double of the beta;
 *   Bhilo    second lowest double of the beta;
 *   Blolo    lowest double of the beta;
 *   Vrehihi  the first dim numbers define the highest doubles of
 *            the real parts of V, the next Householder vector;
 *   Vrelohi  second highest doubles of the real parts of V;
 *   Vrehilo  second lowest doubles of the real parts of V;
 *   Vrelolo  lowest doubles of the real parts of V;
 *   Vimhihi  highest doubles of the imaginary parts of V;
 *   Vimlohi  second highest doubles of the imaginary parts of V;
 *   Vimhilo  second lowest doubles of the imaginary parts of V;
 *   Vimlolo  lowest doubles of the imaginary parts of V;
 *   Wrehihi  the next dim numbers define the highest doubles of
 *            the real parts of the next column in the W matrix;
 *   Wimhihi  the next dim numbers define the highest doubles of
 *            the imaginary parts of the next column in the W matrix;
 *   Wimlohi  the second highest doubles of the imaginary parts of W;
 *   Wimhilo  the second lowest doubles of the imaginary parts of W;
 *   Wimlolo  the lowest doubles of the imaginary parts of W;
 *   WYHrehihi is the dim-by-dim matrix with the highest doubles of the real
 *            parts of W*Y^H;
 *   WYHrelohi are the second highest doubles of the imaginary parts of W*Y^H;
 *   WYHrehilo are the second lowest doubles of the imaginary parts of W*Y^H;
 *   WYHrelolo are the lowest doubles of the imaginary parts of W*Y^H;
 *   WYHimhihi are the highest doubles of the real parts of W*Y^H;
 *   WYHimlohi are the second highest doubles of the real parts of W*Y^H;
 *   WYHimhilo are the second lowest doubles of the imaginary parts of W*Y^H;
 *   WYHimlolo are the lowest doubles of the imaginary parts of W*Y^H.
 *
 * ON RETURN :
 *   Wrehihi  highest doubles of the real parts of the next column of W;
 *   Wrelohi  second highest doubles of the real parts of W;
 *   Wrehilo  second lowest doubles of the real parts of W;
 *   Wrelolo  lowest doubles of the real parts of W;
 *   Wimhihi  highest doubles of the imaginary parts of W;
 *   Wimlohi  second highest doubles of the imaginary parts of W;
 *   Wimhilo  second lowest doubles of the imaginary parts of W;
 *   Wimlolo  lowest doubles of the imaginary parts of W. */

__global__ void dbl4_small_WYT
 ( int nrows, int szt,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo );
/*
 * DESCRIPTION :
 *   Multiplies W with V^T into the matrix WYT, on real data.
 *   Computes Y*W^T, swapping W with Y in the input arguments.
 *
 * ON ENTRY :
 *   nrows    number of rows of all matrices;
 *   szt      number of columns in W and Y,
 *            equals the number of threads in a block;
 *   Whihi    highest doubles of the W matrix in the WY representation;
 *   Wlohi    second highest doubles of W;
 *   Whilo    second lowest doubles of W;
 *   Wlolo    lowest doubles of W;
 *   Vhihi    columns of V are highest doubles of the Householder vectors;
 *   Vlohi    second highest doubles of Householder vectors;
 *   Vhilo    second lowest doubles of Householder vectors;
 *   Vlolo    lowest doubles of Householder vectors;
 *   WYThi    space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTlo    space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt.
 *
 * ON RETURN :
 *   WYThihi  highest doubles of W*Y^T;
 *   WYTlohi  second highest doubles of W*Y^T;
 *   WYThilo  second lowest doubles of W*Y^T;
 *   WYTlolo  lowest doubles of W*Y^T. */

__global__ void cmplx4_small_WYH
 ( int nrows, int szt,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *Yrehihi, double *Yrelohi, double *Yrehilo, double *Yrelolo,
   double *Yimhihi, double *Yimlohi, double *Yimhilo, double *Yimlolo,
   double *WYTrehihi, double *WYTrelohi,
   double *WYTrehilo, double *WYTrelolo,
   double *WYTimhihi, double *WYTimlohi,
   double *WYTimhilo, double *WYTimlolo );
/*
 * DESCRIPTION :
 *   Multiplies W with Y^T into the matrix WYT, on complex data.
 *   Computes Y*W^T, swapping W with Y in the input arguments.
 *
 * ON ENTRY :
 *   nrows    number of rows of all matrices;
 *   szt      number of columns in W and Y,
 *            equals the number of threads in a block;
 *   Wrehihi  highest doubles of the real parts of W;
 *   Wrelohi  second highest doubles of the real parts of W;
 *   Wrehilo  second lowest doubles of the real parts of W;
 *   Wrelolo  lowest doubles of the real parts of W;
 *   Wimhihi  highest doubles of the imaginary parts of W;
 *   Wimlohi  second highest doubles of the imaginary parts of W;
 *   Wimhilo  second lowest doubles of the imaginary parts of W;
 *   Wimlolo  lowest doubles of the imaginary parts of W;
 *   Yrehihi  highest doubles of the real parts of the columns of Y;
 *   Yrelohi  second highest doubles of the real parts of the columns of Y;
 *   Yrehilo  second lowest doubles of the real parts of Y;
 *   Yrelolo  lowest doubles of the real parts of Y;
 *   Yimhihi  highest doubles of the imaginary parts of Y;
 *   Yimlohi  second highest doubles of the imaginary parts of Y;
 *   Yimhilo  second lowest doubles of the imaginary parts of Y;
 *   Yimlolo  lowest doubles of the imaginary parts of Y;
 *   WYTrehihi has space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTrelohi has space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTrehilo has space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTrelolo has space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTimhihi has space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTimlohi has space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTimhilo has space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTimlolo has space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt.
 *
 * ON RETURN :
 *   WYTrehihi are the highest doubles of the real parts of W*Y^H;
 *   WYTrelohi are the second highest doubles of the real parts of W*Y^H;
 *   WYTrelohi are the second lowest doubles of the real parts of W*Y^H;
 *   WYTrelohi are the lowest doubles of the real parts of W*Y^H;
 *   WYTimlohi are the highest doubles of the imaginary parts of W*Y^H;
 *   WYTimlohi are the second highest doubles of the imaginary parts of W*Y^H;
 *   WYTimlohi are the second lowest doubles of the imaginary parts of W*Y^H;
 *   WYTimlohi are the lowest doubles of the imaginary parts of W*Y^H. */

__global__ void dbl4_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihi, double *Qlohi, double *Qhilo, double *Qlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo,
   double *QWYThihi, double *QWYTlohi, double *QWYThilo, double *QWYTlolo );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYT into the matrix QWYT, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYT matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qhihi    highest doubles of Q;
 *   Qlohi    second highest doubles of Q;
 *   Qhilo    second lowest doubles of Q;
 *   Qlolo    lowest doubles of Q;
 *   WYThihi  highest doubles of W*Y^T;
 *   WYTlohi  second highest doubles of W*Y^T;
 *   WYThilo  second lowest doubles of W*Y^T;
 *   WYTlolo  lowest doubles of W*Y^T;
 *   QWYThihi has space for a dim-by-dim matrix;
 *   QWYTlohi has space for a dim-by-dim matrix;
 *   QWYThilo has space for a dim-by-dim matrix;
 *   QWYTlolo has space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYThihi are the highest doubles of Q*WYT;
 *   QWYTlohi are the second highest doubles of Q*WYT;
 *   QWYThilo are the second lowest doubles of Q*WYT;
 *   QWYTlolo are the lowest doubles of Q*WYT. */

__global__ void cmplx4_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihi, double *Qrelohi, double *Qrehilo, double *Qrelolo,
   double *Qimhihi, double *Qimlohi, double *Qimhilo, double *Qimlolo,
   double *WYTrehihi, double *WYTrelohi, double *WYTrehilo, double *WYTrelolo,
   double *WYTimhihi, double *WYTimlohi, double *WYTimhilo, double *WYTimlolo,
   double *QWYTrehihi, double *QWYTrelohi,
   double *QWYTrehilo, double *QWYTrelolo,
   double *QWYTimhihi, double *QWYTimlohi,
   double *QWYTimhilo, double *QWYTimlolo );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYT into the matrix QWYT, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYT matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qrehihi  highest doubles of the real parts of Q;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   WYTrehihi are the highest doubles of the real parts of W*Y^T;
 *   WYTrelohi are the second highest doubles of the real parts of W*Y^T;
 *   WYTrehilo are the second lowest doubles of the real parts of W*Y^T;
 *   WYTrelolo are the lowest doubles of the real parts of W*Y^T;
 *   WYTimhihi are the highest doubles of the imaginary parts of W*Y^T;
 *   WYTimlohi are the second highest doubles of the imaginary parts of W*Y^T;
 *   WYTimhilo are the second lowest doubles of the imaginary parts of W*Y^T;
 *   WYTimlolo are the lowest doubles of the imaginary parts of W*Y^T;
 *   QWYTrehihi has space for a dim-by-dim matrix;
 *   QWYTrelohi has space for a dim-by-dim matrix;
 *   QWYTrehilo has space for a dim-by-dim matrix;
 *   QWYTrelolo has space for a dim-by-dim matrix;
 *   QWYTimhihi has space for a dim-by-dim matrix;
 *   QWYTimlohi has space for a dim-by-dim matrix;
 *   QWYTimhilo has space for a dim-by-dim matrix;
 *   QWYTimlolo has space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYTrehihi are the highest doubles of the real parts of Q*W*Y^T;
 *   QWYTrelohi are the second highest doubles of the real parts of Q*W*Y^T;
 *   QWYTrehilo are the second lowest doubles of the real parts of Q*W*Y^T;
 *   QWYTrelolo are the lowest doubles of the real parts of Q*W*Y^T;
 *   QWYTimhihi are the highest doubles of the imaginary parts of Q*W*Y^T;
 *   QWYTimlohi are the second highest doubles of the imaginary parts
 *              of Q*W*Y^T;
 *   QWYTimhilo are the second lowest doubles of the imaginary parts
 *              of Q*W*Y^T;
 *   QWYTimlolo are the lowest doubles of the imaginary parts of Q*W*Y^T. */

__global__ void dbl4_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWThihi, double *YWTlohi, double *YWThilo, double *YWTlolo,
   double *Chihi, double *Clohi, double *Chilo, double *Clolo,
   double *YWTChihi, double *YWTClohi, double *YWTChilo, double *YWTClolo );
/*
 * DESCRIPTION :
 *   Multiplies YWT with C into the matrix YWTC, on real data.
 *
 * ON ENTRY :
 *   nrows    total number of rows, nrows = rowdim + rowoff;
 *   ncols    total number of colums, ncols = coldim + coloff;
 *   rowdim   number of rows minus the row offset;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in C;
 *   coloff   offset for the column index in C;
 *   Chihi    highest doubles of C;
 *   Clohi    second highest doubles of C;
 *   Chilo    second lowest doubles of C;
 *   Clolo    lowest doubles of C;
 *   YWThihi  highest doubles of Y*W^T;
 *   YWTlohi  second highest doubles of Y*W^T;
 *   YWThilo  second lowest doubles of Y*W^T;
 *   YWTlolo  lowest doubles of Y*W^T;
 *   YWTChihi has space for a rowdim-by-coldim matrix;
 *   YWTClohi has space for a rowdim-by-coldim matrix;
 *   YWTChilo has space for a rowdim-by-coldim matrix;
 *   YWTClolo has space for a rowdim-by-coldim matrix.
 *
 * ON RETURN :
 *   YWTChihi are the highest doubles of the product of YWT with C;
 *   YWTClohi are the second highest doubles of the product of YWT with C;
 *   YWTChilo are the second lowest doubles of the product of YWT with C;
 *   YWTClolo are the lowest doubles of the product of YWT with C. */

__global__ void cmplx4_small_YWHC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWTrehihi, double *YWTrelohi, double *YWTrehilo, double *YWTrelolo,
   double *YWTimhihi, double *YWTimlohi, double *YWTimhilo, double *YWTimlolo,
   double *Crehihi, double *Crelohi, double *Crehilo, double *Crelolo,
   double *Cimhihi, double *Cimlohi, double *Cimhilo, double *Cimlolo,
   double *YWTCrehihi, double *YWTCrelohi,
   double *YWTCrehilo, double *YWTCrelolo,
   double *YWTCimhihi, double *YWTCimlohi,
   double *YWTCimhilo, double *YWTCimlolo );
/*
 * DESCRIPTION :
 *   Multiplies YWT with C into the matrix YWTC, on complex data.
 *
 * ON ENTRY :
 *   nrows    total number of rows, nrows = rowdim + rowoff;
 *   ncols    total number of colums, ncols = coldim + coloff;
 *   rowdim   number of rows minus the row offset;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in C;
 *   coloff   offset for the column index in C;
 *   Crehihi  highest doubles of the real parts of C;
 *   Crelohi  second highest doubles of the real parts of C;
 *   Crehilo  second lowest doubles of the real parts of C;
 *   Crelolo  lowest doubles of the real parts of C;
 *   Cimhihi  highest doubles of the imaginary parts of C;
 *   Cimlohi  second highest doubles of the imaginary parts of C;
 *   Cimhilo  second lowest doubles of the imaginary parts of C;
 *   Cimlolo  lowest doubles of the imaginary parts of C;
 *   YWTrehihi has the highest doubles of the real parts of Y*W^T;
 *   YWTrelohi has the second highest doubles of the real parts of Y*W^T;
 *   YWTrehilo has the second lowest doubles of the real parts of Y*W^T;
 *   YWTrelolo has the lowest doubles of the real parts of Y*W^T;
 *   YWTimhihi has the highest doubles of the imaginary parts of Y*W^T;
 *   YWTimlohi has the second highest doubles of the imaginary parts of Y*W^T;
 *   YWTimhilo has the second lowest doubles of the imaginary parts of Y*W^T;
 *   YWTimlolo has the lowest doubles of the imaginary parts of Y*W^T;
 *   YWTCrehihi has space for a rowdim-by-coldim matrix;
 *   YWTCrelohi has space for a rowdim-by-coldim matrix;
 *   YWTCrehilo has space for a rowdim-by-coldim matrix;
 *   YWTCrelolo has space for a rowdim-by-coldim matrix;
 *   YWTCimhihi has space for a rowdim-by-coldim matrix.
 *   YWTCimlohi has space for a rowdim-by-coldim matrix;
 *   YWTCimhilo has space for a rowdim-by-coldim matrix;
 *   YWTCimlolo has space for a rowdim-by-coldim matrix.
 *
 * ON RETURN :
 *   YWTCrehihi are the highest doubles of the real parts of YWT*C;
 *   YWTCrelohi are the second highest doubles of the real parts of YWT*C;
 *   YWTCrehilo are the second lowest doubles of the real parts of YWT*C;
 *   YWTCrelolo are the lowest doubles of the real parts of YWT*C;
 *   YWTCimhihi are the highest doubles of the imaginary parts of YWT*C;
 *   YWTCimlohi are the second highest doubles of the imaginary parts
 *            of YWT*C;
 *   YWTCimhilo are the second lowest doubles of the imaginary parts
 *            of YWT*C;
 *   YWTCimlolo are the lowest doubles of the imaginary parts of YWT*C. */

__global__ void dbl4_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihi, double *Qlohi, double *Qhilo, double *Qlolo,
   double *QWYThihi, double *QWYTlohi, double *QWYThilo, double *QWYTlolo );
/*
 * DESCRIPTION :
 *   Updates Q by adding QWYT, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qhihi    highest doubles of the current orthogonal matrix;
 *   Qlohi    second highest doubles of the current orthogonal matrix;
 *   Qhilo    second lowest doubles of the current orthogonal matrix;
 *   Qlolo    lowest doubles of the current orthogonal matrix;
 *   QWYThihi has space for a dim-by-dim matrix;
 *   QWYTlohi has space for a dim-by-dim matrix;
 *   QWYThilo has space for a dim-by-dim matrix;
 *   QWYTlolo has space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of Q + QWYT;
 *   Qlohi    second highest doubles of Q + QWYT;
 *   Qhilo    second lowest doubles of Q + QWYT;
 *   Qlolo    lowest doubles of Q + QWYT. */

__global__ void cmplx4_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihi, double *Qrelohi, double *Qrehilo, double *Qrelolo,
   double *Qimhihi, double *Qimlohi, double *Qimhilo, double *Qimlolo,
   double *QWYTrehihi, double *QWYTrelohi,
   double *QWYTrehilo, double *QWYTrelolo,
   double *QWYTimhihi, double *QWYTimlohi,
   double *QWYTimhilo, double *QWYTimlolo );
/*
 * DESCRIPTION :
 *   Updates Q by adding QWYT, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   rowdim   dimension minus the column offset;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qrehihi  highest doubles of the real parts of Q;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   QWYTrehihi has space for a dim-by-dim matrix;
 *   QWYTrelohi has space for a dim-by-dim matrix;
 *   QWYTrehilo has space for a dim-by-dim matrix;
 *   QWYTrelolo has space for a dim-by-dim matrix;
 *   QWYTimhihi has space for a dim-by-dim matrix;
 *   QWYTimlohi has space for a dim-by-dim matrix;
 *   QWYTimhilo has space for a dim-by-dim matrix;
 *   QWYTimlolo has space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of Q + QWYT;
 *   Qrelohi  second highest doubles of the real parts of Q + QWYT;
 *   Qrehilo  second lowest doubles of the real parts of Q + QWYT;
 *   Qrelolo  lowest doubles of the real parts of Q + QWYT;
 *   Qimhihi  highest doubles of the imaginary parts Q + QWYT;
 *   Qimlohi  second highest doubles of the imaginary parts Q + QWYT;
 *   Qimhilo  second lowest doubles of the imaginary parts Q + QWYT;
 *   Qimlolo  lowest doubles of the imaginary parts Q + QWYT. */

__global__ void dbl4_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *YWTChihi, double *YWTClohi, double *YWTChilo, double *YWTClolo );
/*
 * DESCRIPTION :
 *   Updates R by adding YWTC, on real data.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in R;
 *   coloff   offset for the column index in R;
 *   Rhihi    highest doubles of R;
 *   Rlohi    second highest doubles of R;
 *   Rhilo    second lowest doubles of R;
 *   Rlolo    lowest doubles of R;
 *   YWTChihi are the highest doubles of YWT*C;
 *   YWTClohi are the second highest doubles of YWT*C;
 *   YWTChilo are the second lowest doubles of YWT*C;
 *   YWTClolo are the lowest doubles of YWT*C.
 *
 * ON RETURN :
 *   Rhihi    highest doubles of R + YWTC;
 *   Rlohi    second highest doubles of R + YWTC;
 *   Rhilo    second lowest doubles of R + YWTC;
 *   Rlolo    lowest doubles of R + YWTC. */

__global__ void cmplx4_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *YWTCrehihi, double *YWTCrelohi,
   double *YWTCrehilo, double *YWTCrelolo,
   double *YWTCimhihi, double *YWTCimlohi,
   double *YWTCimhilo, double *YWTCimlolo );
/*
 * DESCRIPTION :
 *   Updates R by adding YWTC, on complex data.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in R;
 *   coloff   offset for the column index in R;
 *   Rrehi    highest doubles of the real parts of R;
 *   Rrehi    second highest doubles of the real parts of R;
 *   Rrelo    second lowest doubles of the real parts of R;
 *   Rrelo    lowest doubles of the real parts of R;
 *   Rimhi    highest doubles of the imaginary parts of R;
 *   Rimhi    second highest doubles of the imaginary parts of R;
 *   Rimlo    second lowest doubles of the imaginary parts of R;
 *   Rimlo    lowest doubles of the imaginary parts of R;
 *   YWTCrehihi are the highest doubles of the real parts of YWT*C;
 *   YWTCrelohi are the second highest doubles of the real parts of YWT*C;
 *   YWTCrehilo are the second lowest doubles of the real parts of YWT*C;
 *   YWTCrelolo are the lowest doubles of the real parts of YWT*C;
 *   YWTCimhihi are the highest doubles of the imaginary parts of YWT*C;
 *   YWTCimlohi are the second highest doubles of the imaginary parts
 *            of YWT*C;
 *   YWTCimhilo are the second lowest doubles of the imaginary parts
 *            of YWT*C;
 *   YWTCimlolo are the lowest doubles of the imaginary parts of YWT*C.
 *
 * ON RETURN :
 *   Rrehihi  highest doubles of the real parts of R + YWTC;
 *   Rrelohi  second highest doubles of the real parts of R + YWTC;
 *   Rrehilo  second lowest doubles of the real parts of R + YWTC;
 *   Rrelolo  lowest doubles of the real parts of R + YWTC;
 *   Rimhihi  highest doubles of the imaginary parts R + YWTC;
 *   Rimlohi  second highest doubles of the imaginary parts R + YWTC;
 *   Rimhilo  second lowest doubles of the imaginary parts R + YWTC;
 *   Rimlolo  lowest doubles of the imaginary parts R + YWTC. */

void GPU_dbl4_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *vhihi_h, double *vlohi_h, double *vhilo_h, double *vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the Householder vector for small
 *   enough matrices for the vector to fit entirely in shared memory.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   nbt      number of tiles, szt*nbt = ncols;
 *   colidx   global index of the current column;
 *   nrows1   number of threads in the block equals the number
 *            of elements computed in the Householder vector;
 *   L        local index of the column in the current tile;
 *   Ahihi_h  highest doubles of the matrix on the host;
 *   Alohi_h  second highest doubles of the matrix on the host;
 *   Ahilo_h  second lowest doubles of the matrix on the host;
 *   Alolo_h  lowest doubles of the matrix on the host;
 *   Ahihi_d  highest doubles of the matrix on the device;
 *   Alohi_d  second highest doubles of the matrix on the device;
 *   Ahilo_d  secdon lowest doubles of the matrix on the device;
 *   Alolo_d  lowest doubles of the matrix on the device;
 *   vhihi_h  space for the current Householder vector;
 *   vlohi_h  space for the current Householder vector;
 *   vhilo_h  space for the current Householder vector;
 *   vlolo_h  space for the current Householder vector;
 *   Vhihi_d  space for the Householder vectors on the device;
 *   Vlohi_d  space for the Householder vectors on the device;
 *   Vhilo_d  space for the Householder vectors on the device;
 *   Vlolo_d  space for the Householder vectors on the device;
 *   betahihi_h has space for the highest doubles of the betas, if verbose;
 *   betalohi_h has space for the second highest doubles of the betas,
 *            if verbose;
 *   betahilo_h has space for the second lowest doubles of the betas,
 *            if verbose;
 *   betalolo_h has space for the lowest doubles of the betas, if verbose;
 *   betahihi_d has space on the device for the highest doubles of the betas;
 *   betalohi_d has space on the device for the second highest doubles
 *            of the betas;
 *   betahilo_d has space on the device for the second lowest doubles
 *            of the betas;
 *   betalolo_d has space on the device for the lowest doubles of the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   sqrtfun  current number of calls to sqrt();
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   vhihi_h  highest doubles of v,
 *            the next Householder vector on the host, if verbose;
 *   vlohi_h  second highest doubles of v;
 *   vhilo_h  second lowest doubles of v;
 *   vlolo_h  lowest doubles of v;
 *   Vhihi_d  highest doubles of the next computed Householder vector;
 *   Vlohi_d  second highest doubles of the next computed Householder vector;
 *   Vhilo_d  second lowest doubles of the next computed Householder vector;
 *   Vlolo_d  lowest doubles of the next computed Householder vector;
 *   betahihi_h has the updated vector of the highest doubles of the betas,
 *            if verbose;
 *   betalohi_h has the updated vector of the second highest doubles
 *            of the betas, if verbose;
 *   betahilo_h has the updated vector of the second lowest doubles
 *            of the betas, if verbose;
 *   betalolo_h has the updated vector of the lowest doubles of the betas,
 *            if verbose;
 *   betahihi_d has the highest double of the next beta constant;
 *   betalohi_d has the second highest double of the next beta constant;
 *   betahilo_d has the second lowest double of the next beta constant;
 *   betalolo_d has the lowest double of the next beta constant;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions;
 *   sqrtfun  accumulated number of calls to sqrt(). */

void GPU_cmplx4_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *vrehihi_h, double *vrelohi_h, double *vrehilo_h, double *vrelolo_h,
   double *vimhihi_h, double *vimlohi_h, double *vimhilo_h, double *vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the Householder vector for small
 *   enough matrices for the vector to fit entirely in shared memory.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   nbt      number of tiles, szt*nbt = ncols;
 *   colidx   global index of the current column;
 *   nrows1   number of threads in the block equals the number
 *            of elements computed in the Householder vector;
 *   L        local index of the column in the current tile;
 *   Arehihi_h are the highest doubles of the real parts of A on the host;
 *   Arelohi_h are the second highest doubles of the real parts of A
 *            on the host;
 *   Arehilo_h are the second lowest doubles of the real parts of A
 *            on the host;
 *   Arelolo_h are the lowest doubles of the real parts of A on the host;
 *   Aimhihi_h are the highest doubles of the imaginary parts of A
 *            on the host;
 *   Aimlohi_h are the second highest doubles of the imaginary parts of A
 *            on the host;
 *   Aimhilo_h are the second lowest doubles of the imaginary parts of A
 *            on the host;
 *   Aimlolo_h are the lowest doubles of the imaginary parts of A on the host;
 *   Arehihi_d are the highest doubles of the real parts of A on the device;
 *   Arelohi_d are the second highest doubles of the real parts of A
 *            on the device;
 *   Arehilo_d are the second lowest doubles of the real parts of A
 *            on the device;
 *   Arelolo_d are the lowest doubles of the real parts of A on the device;
 *   Aimhihi_d are the highest doubles of the imaginary parts of A
 *            on the device;
 *   Aimlohi_d are the second highest doubles of the imaginary parts of A,
 *            on the device;
 *   Aimhilo_d are the second lowest doubles of the imaginary parts of A,
 *            on the device;
 *   Aimlolo_d are the lowest doubles of the imaginary parts of A
 *            on the device;
 *   vrehihi_h has space for the highest doubles of the real parts
 *            of the current Householder vector;
 *   vrelohi_h has space for the second highest doubles of the real parts
 *            of the current Householder vector;
 *   vrehilo_h has space for the second lowest doubles of the real parts
 *            of the current Householder vector;
 *   vrelolo_h has space for the lowest doubles of the real parts
 *            of the current Householder vector;
 *   vimhihi_h has space for the highest doubles of the imaginary parts
 *            of the current Householder vector;
 *   vimlohi_h has space for the second highest doubles of the imaginary
 *            parts of the current Householder vector;
 *   vimhilo_h has space for the second lowest doubles of the imaginary parts
 *            of the current Householder vector;
 *   vimlolo_h has space for the lowest doubles of the imaginary parts
 *            of the current Householder vector;
 *   Vrehihi_d has space for the highest doubles of the real parts
 *            of the Householder vectors on the device;
 *   Vrelohi_d has space for the second highest doubles of the real parts
 *            of the Householder vectors on the device;
 *   Vrehilo_d has space for the second lowest doubles of the real parts
 *            of the Householder vectors on the device;
 *   Vrelolo_d has space for the lowest doubles of the real parts
 *            of the Householder vectors on the device;
 *   Vimhihi_d has space for the highest doubles of the imaginary parts
 *            of the Householder vectors on the device;
 *   Vimlohi_d has space for the second highest doubles of the imaginary
 *            parts of the Householder vectors on the device;
 *   Vimhilo_d has space for the second lowest doubles of the imaginary parts
 *            of the Householder vectors on the device;
 *   Vimlolo_d has space for the lowest doubles of the imaginary parts
 *            of the Householder vectors on the device;
 *   betahihi_h has space for the highest doubles of betas if verbose;
 *   betalohi_h has space for the second highest doubles of betas if verbose;
 *   betahilo_h has space for the second lowest doubles of betas if verbose;
 *   betalolo_h has space for the lowest doubles of betas if verbose;
 *   betahihi_d has space on the device for the highest doubles of the betas;
 *   betalohi_d has space on the device for the second highest doubles
 *            of the betas;
 *   betahilo_d has space on the device for the second lowest doubles
 *            of the betas;
 *   betalolo_d has space on the device for the lowest doubles of the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   sqrtfun  current number of calls to sqrt();
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   vrehihi_h are the highest doubles of the real parts of v,
 *            the next Householder vector on the host, if verbose;
 *   vrelohi_h are the second highest doubles of the real parts of v;
 *   vrehilo_h are the second lowest doubles of the real parts of v;
 *   vrelolo_h are the lowest doubles of the real parts of v;
 *   Vrehihi_d are the highest doubles of the real parts 
 *            of the next computed Householder vector V;
 *   Vrelohi_d are the second highest doubles of the real parts of V;
 *   Vrehilo_d are the second lowest doubles of the real parts of V;
 *   Vrelolo_d are the lowest doubles of the real parts of V;
 *   Vimhihi_d are the highest doubles of the imaginary parts of V;
 *   Vimlohi_d are the second highest doubles of the imaginary parts of V;
 *   Vimhilo_d are the second lowest doubles of the imaginary parts of V;
 *   Vimlolo_d are the lowest doubles of the imaginary parts of V;
 *   betahi_h has the highest doubles of the updated vector
 *            of betas, if verbose;
 *   betahi_h has the second highest doubles of the updated vector
 *            of betas, if verbose;
 *   betalo_h has the second lowest doubles of the updated vector
 *            of betas, if verbose;
 *   betalo_h has the lowest doubles of the updated vector
 *            of betas, if verbose;
 *   betahi_d has the highest doubles of the next beta constant;
 *   betahi_d has the second highest doubles of the next beta constant;
 *   betalo_d has the second lowest doubles of the next beta constant;
 *   betalo_d has the lowest doubles of the next beta constant;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions;
 *   sqrtfun  accumulated number of calls to sqrt(). */

void GPU_dbl4_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *vhihi_h, double *vlohi_h, double *vhilo_h, double *vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *sumshihi_h, double *sumslohi_h,
   double *sumshilo_h, double *sumslolo_h,
   double *sumshihi_d, double *sumslohi_d,
   double *sumshilo_d, double *sumslolo_d,
   double *sigmahihi_h, double *sigmalohi_h, 
   double *sigmahilo_h, double *sigmalolo_h, 
   double *sigmahihi_d, double *sigmalohi_d,
   double *sigmahilo_d, double *sigmalolo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernels to compute the Householder vector for matrices
 *   of any size, with multiple blocks of threads, on real data.
 *
 * REQUIRED : nrows1 > szt, for nrows1 <= szt, call GPU_dbl4_small_house.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile, equals the number of threads in a block;
 *   nbt      number of tiles, szt*nbt = ncols;
 *   colidx   global index of the current column;
 *   nrows1   number of threads in the block equals the number
 *            of elements computed in the Householder vector;
 *   L        local index of the column in the current tile;
 *   Ahihi_h  highest doubles of the matrix on the host;
 *   Alohi_h  second highest doubles of the matrix on the host;
 *   Ahilo_h  second lowest doubles of the matrix on the host;
 *   Alolo_h  lowest doubles of the matrix on the host;
 *   Ahihi_d  highest doubles of the matrix on the device;
 *   Alohi_d  second highest doubles of the matrix on the device;
 *   Ahilo_d  secdon lowest doubles of the matrix on the device;
 *   Alolo_d  lowest doubles of the matrix on the device;
 *   vhihi_h  space for the current Householder vector, on the host;
 *   vlohi_h  space for the current Householder vector, on the host;
 *   vhilo_h  space for the current Householder vector, on the host;
 *   vlolo_h  space for the current Householder vector, on the host;
 *   Vhihi_d  space for the Householder vectors, on the device;
 *   Vlohi_d  space for the Householder vectors, on the device;
 *   Vhilo_d  space for the Householder vectors, on the device;
 *   Vlolo_d  space for the Householder vectors, on the device;
 *   betahihi_h has space for the highest doubles of the betas, if verbose;
 *   betalohi_h has space for the second highest doubles of the betas,
 *            if verbose;
 *   betahilo_h has space for the second lowest doubles of the betas,
 *            if verbose;
 *   betalolo_h has space for the lowest doubles of the betas if verbose;
 *   betahihi_d has space on the device for the highest doubles of the betas;
 *   betalohi_d has space on the device for the second highest doubles
 *            of the betas;
 *   betahilo_d has space on the device for the second lowest doubles
 *            of the betas;
 *   betalolo_d has space on the device for the lowest doubles of the betas;
 *   sumshihi_h has space for the highest doubles of sums, if verbose;
 *   sumslohi_h has space for the second highest doubles of sums, if verbose;
 *   sumshilo_h has space for the second lowest doubles of sums, if verbose;
 *   sumslolo_h has space for the lowest doubles of sums, if verbose;
 *   sumshihi_d has space for the highest doubles of sums, on the device;
 *   sumslohi_d has space for the second highest doubles of sums,
 *            on the device;
 *   sumshilo_d has space for the second lowest doubles of sums,
 *            on the device;
 *   sumslolo_d has space for the lowest doubles of sums, on the device;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   sqrtfun  current number of calls to sqrt();
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   vhihi_h  highest doubles of the next Householder vector v
 *            on the host, if verbose;
 *   vlohi_h  second highest doubles of v;
 *   vhilo_h  second lowest doubles of v;
 *   vlolo_h  lowest doubles of v;
 *   Vhihi_d  highest doubles of the next computed Householder vector V;
 *   Vlohi_d  second highest doubles of V;
 *   Vhilo_d  second lowest doubles of V;
 *   Vlolo_d  lowest doubles of V;
 *   betahihi_h has the updated vector of the highest doubles
 *            of the betas, if verbose;
 *   betalohi_h has the updated vector of the second highest doubles
 *            of the betas, if verbose;
 *   betahilo_h has the updated vector of the second lowest doubles
 *            of the betas, if verbose;
 *   betalolo_h has the updated vector of the lowest doubles
 *            of the betas, if verbose;
 *   betahihi_d has the highest double of the next beta constant;
 *   betalohi_d has the second highest double of the next beta constant;
 *   betahilo_d has the second lowest double of the next beta constant;
 *   betalolo_d has the lowest double of the next beta constant;
 *   sumshihi_h are the highest doubles of the sums, if verbose;
 *   sumslohi_h are the second highest doubles of the sums, if verbose;
 *   sumshilo_h are the second lowest doubles of the sums, if verbose;
 *   sumslolo_h are the lowest doubles of the sums, if verbose;
 *   sumshihi_h are the highest doubles of the sums, on the device;
 *   sumslohi_h are the second highest doubles of the sums, on the device;
 *   sumshilo_h are the second lowest doubles of the sums, on the device;
 *   sumslolo_h are the lowest doubles of the sums, on the device;
 *   sigmahihi_h is the highest double of sigma, on the host;
 *   sigmalohi_h is the second highest double of sigma, on the host;
 *   sigmahilo_h is the second lowest double of sigma, on the host;
 *   sigmalolo_h is the lowest double of sigma, on the host;
 *   sigmahihi_d is the highest double of sigma, on the device;
 *   sigmalohi_d is the second highest double of sigma, on the device;
 *   sigmahilo_d is the second lowest double of sigma, on the device;
 *   sigmalolo_d is the lowest double of sigma, on the device;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions;
 *   sqrtfun  accumulated number of calls to sqrt(). */

void GPU_cmplx4_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *vrehihi_h, double *vrelohi_h, double *vrehilo_h, double *vrelolo_h,
   double *vimhihi_h, double *vimlohi_h, double *vimhilo_h, double *vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *sumshihi_h, double *sumslohi_h,
   double *sumshilo_h, double *sumslolo_h,
   double *sumshihi_d, double *sumslohi_d,
   double *sumshilo_d, double *sumslolo_d,
   double *sigmahihi_h, double *sigmalohi_h,
   double *sigmahilo_h, double *sigmalolo_h,
   double *sigmahihi_d, double *sigmalohi_d,
   double *sigmahilo_d, double *sigmalolo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernels to compute the Householder vector for matrices
 *   of any size, with multiple blocks of threads, on real data.
 *
 * REQUIRED : nrows1 > szt, for nrows1 <= szt, call GPU_dbl4_small_house.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile, equals the number of threads in a block;
 *   nbt      number of tiles, szt*nbt = ncols;
 *   colidx   global index of the current column;
 *   nrows1   number of threads in the block equals the number
 *            of elements computed in the Householder vector;
 *   L        local index of the column in the current tile;
 *   Arehihi_h are the highest doubles of the real parts of A on the host;
 *   Arelohi_h are the second highest doubles of the real parts of A
 *            on the host;
 *   Arehilo_h are the second lowest doubles of the real parts of A
 *            on the host;
 *   Arelolo_h are the lowest doubles of the real parts of A on the host;
 *   Aimhihi_h are the highest doubles of the imaginary parts of A
 *            on the host;
 *   Aimlohi_h are the second highest doubles of the imaginary parts of A
 *            on the host;
 *   Aimhilo_h are the second lowest doubles of the imaginary parts of A
 *            on the host;
 *   Aimlolo_h are the lowest doubles of the imaginary parts of A on the host;
 *   Arehihi_d are the highest doubles of the real parts of A on the device;
 *   Arelohi_d are the second highest doubles of the real parts of A
 *            on the device;
 *   Arehilo_d are the second lowest doubles of the real parts of A
 *            on the device;
 *   Arelolo_d are the lowest doubles of the real parts of A on the device;
 *   Aimhihi_d are the highest doubles of the imaginary parts of A
 *            on the device;
 *   Aimlohi_d are the second highest doubles of the imaginary parts of A,
 *            on the device;
 *   Aimhilo_d are the second lowest doubles of the imaginary parts of A,
 *            on the device;
 *   Aimlolo_d are the lowest doubles of the imaginary parts of A
 *            on the device;
 *   vrehihi_h has space for the current Householder vector;
 *   vrelohi_h has space for the current Householder vector;
 *   vrehilo_h has space for the current Householder vector;
 *   vrelolo_h has space for the current Householder vector;
 *   vimhihi_h has space for the current Householder vector;
 *   vimlohi_h has space for the current Householder vector;
 *   vimhilo_h has space for the current Householder vector;
 *   vimlolo_h has space for the current Householder vector;
 *   Vrehihi_d has space for the Householder vectors on the device;
 *   Vrelohi_d has space for the Householder vectors on the device;
 *   Vrehilo_d has space for the Householder vectors on the device;
 *   Vrelolo_d has space for the Householder vectors on the device;
 *   Vimhihi_d has space for the Householder vectors on the device;
 *   Vimlohi_d has space for the Householder vectors on the device;
 *   Vimhilo_d has space for the Householder vectors on the device;
 *   Vimlolo_d has space for the Householder vectors on the device;
 *   betahihi_h has space for the highest doubles of the betas if verbose;
 *   betalohi_h has space for the second highest doubles of the betas,
 *            if verbose;
 *   betahilo_h has space for the second lowest doubles of the betas,
 *            if verbose;
 *   betalolo_h has space for the lowest doubles of the betas if verbose;
 *   betahihi_d has space on the device for the highest doubles of the betas;
 *   betalohi_d has space on the device for the second highest doubles
 *            of the betas;
 *   betahilo_d has space on the device for the second lowest doubles
 *            of the betas;
 *   betalolo_d has space on the device for the lowest doubles of the betas;
 *   sumshihi_h has space for the highest doubles of sums, if verbose;
 *   sumslohi_h has space for the second highest doubles of sums, if verbose;
 *   sumshilo_h has space for the second lowest doubles of sums, if verbose;
 *   sumslolo_h has space for the lowest doubles of sums, if verbose;
 *   sumshihi_d has space for the highest doubles of sums, on the device;
 *   sumslohi_d has space for the second highest doubles of sums,
 *            on the device;
 *   sumshilo_d has space for the second lowest doubles of sums,
 *            on the device;
 *   sumslolo_d has space for the lowest doubles of sums, on the device;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   sqrtfun  current number of calls to sqrt();
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   vrehihi_h are the highest doubles of the real parts of v,
 *            the next Householder vector on the host, if verbose;
 *   vrelohi_h are the second highest doubles of the real parts of v;
 *   vrehilo_h are the second lowest doubles of the real parts of v;
 *   vrelolo_h are the lowest doubles of the real parts of v;
 *   vimhihi_h are the highest doubles of imaginary parts of ;
 *   vimlohi_h are the second highest doubles of imaginary parts of v;
 *   vimhilo_h are the second lowest doubles of imaginary parts of v;
 *   vimlolo_h are the lowest doubles of imaginary parts of v;
 *   Vrehihi_d are the highest doubles of the real parts of V,
 *            the next Householder vector, on the device;
 *   Vrehilo_d are the second highest doubles of the real parts of V;
 *   Vrehilo_d are the second lowest doubles of the real parts of V;
 *   Vrelolo_d are the lowest doubles of the real parts of V;
 *   Vimhihi_d are the highest doubles of the imaginary parts of V;
 *   Vimlohi_d are the second highest doubles of the imaginary parts of V;
 *   Vimhilo_d are the second lowest doubles of the imaginary parts of V;
 *   Vimlolo_d are the lowest doubles of the imaginary parts of V;
 *   betahihi_h has the updated vector of the highest doubles
 *            of the betas if verbose;
 *   betalohi_h has the updated vector of the second highest doubles
 *            of the betas if verbose;
 *   betahilo_h has the updated vector of the second lowest doubles
 *            of the betas if verbose;
 *   betalolo_h has the updated vector of the lowest doubles
 *            of the betas if verbose;
 *   betahihi_d has the highest double of the next beta constant;
 *   betalohi_d has the second highest double of the next beta constant;
 *   betahilo_d has the second lowest double of the next beta constant;
 *   betalolo_d has the lowest double of the next beta constant;
 *   sumshihi_h are the highest doubles of the sums, if verbose;
 *   sumslohi_h are the second highest doubles of the sums, if verbose;
 *   sumshilo_h are the second lowest doubles of the sums, if verbose;
 *   sumslolo_h are the lowest doubles of the sums, if verbose;
 *   sumshihi_h are the highest doubles of the sums, on the device;
 *   sumslohi_h are the second highest doubles of the sums, on the device;
 *   sumshilo_h are the second lowest doubles of the sums, on the device;
 *   sumslolo_h are the lowest doubles of the sums, on the device;
 *   sigmahihi_h is the highest double of sigma, on the host;
 *   sigmalohi_h is the second highest double of sigma, on the host;
 *   sigmahilo_h is the second lowest double of sigma, on the host;
 *   sigmalolo_h is the lowest double of sigma, on the host;
 *   sigmahihi_d is the highest double of sigma, on the device;
 *   sigmalohi_d is the second highest double of sigma, on the device;
 *   sigmahilo_d is the second lowest double of sigma, on the device;
 *   sigmalolo_d is the lowest double of sigma, on the device;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions;
 *   sqrtfun  accumulated number of calls to sqrt(). */

void GPU_dbl4_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update one tile.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the reduced matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   colidx   global index of the current column;
 *   k        index of the current tile;
 *   L        local index of the column in the current tile;
 *   Ahihi_h  highest doubles of the matrix on the host;
 *   Alohi_h  second highest doubles of the matrix on the host;
 *   Ahilo_h  second lowest doubles of the matrix on the host;
 *   Alolo_h  lowest doubles of the matrix on the host;
 *   Ahihi_d  highest doubles of the matrix on the device;
 *   Alohi_d  second highest doubles of the matrix on the device;
 *   Ahilo_d  secdon lowest doubles of the matrix on the device;
 *   Alolo_d  lowest doubles of the matrix on the device;
 *   Vhihi_d  space allocated for the Householder vectors on the device;
 *   Vlohi_d  space allocated for the Householder vectors on the device;
 *   Vhilo_d  space allocated for the Householder vectors on the device;
 *   Vlolo_d  space allocated for the Householder vectors on the device;
 *   betahihi_h has space allocated for the betas if verbose;
 *   betalohi_h has space allocated for the betas if verbose;
 *   betahilo_h has space allocated for the betas if verbose;
 *   betalolo_h has space allocated for the betas if verbose;
 *   betahihi_d has space allocated on the device for the betas;
 *   betalohi_d has space allocated on the device for the betas;
 *   betahilo_d has space allocated on the device for the betas;
 *   betalolo_d has space allocated on the device for the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Ahihi_d  highest doubles of the reduced matrix on the device;
 *   Alohi_d  second highest doubles of the reduced matrix on the device;
 *   Ahilo_d  second lowest doubles of the reduced matrix on the device;
 *   Alolo_d  lowest doubles of the reduced matrix on the device;
 *   Ahihi_h  highest doubles of the reduced matrix on the host,
 *            if verbose;
 *   Alohi_h  second highest doubles of the reduced matrix on the host,
 *            if verbose;
 *   Ahilo_h  second lowest doubles of the reduced matrix on the host,
 *            if verbose;
 *   Alolo_h  lowest doubles of the reduced matrix on the host,
 *            if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx4_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update one tile.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the reduced matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   colidx   global index of the current column;
 *   k        index of the current tile;
 *   L        local index of the column in the current tile;
 *   Arehihi_h are the highest doubles of the real parts of A on the host;
 *   Arelohi_h are the second highest doubles of the real parts of A
 *            on the host;
 *   Arehilo_h are the second lowest doubles of the real parts of A
 *            on the host;
 *   Arelolo_h are the lowest doubles of the real parts of A on the host;
 *   Aimhihi_h are the highest doubles of the imaginary parts of A
 *            on the host;
 *   Aimlohi_h are the second highest doubles of the imaginary parts of A
 *            on the host;
 *   Aimhilo_h are the second lowest doubles of the imaginary parts of A
 *            on the host;
 *   Aimlolo_h are the lowest doubles of the imaginary parts of A on the host;
 *   Arehihi_d are the highest doubles of the real parts of A on the device;
 *   Arelohi_d are the second highest doubles of the real parts of A
 *            on the device;
 *   Arehilo_d are the second lowest doubles of the real parts of A
 *            on the device;
 *   Arelolo_d are the lowest doubles of the real parts of A on the device;
 *   Aimhihi_d are the highest doubles of the imaginary parts of A
 *            on the device;
 *   Aimlohi_d are the second highest doubles of the imaginary parts of A,
 *            on the device;
 *   Aimhilo_d are the second lowest doubles of the imaginary parts of A,
 *            on the device;
 *   Aimlolo_d are the lowest doubles of the imaginary parts of A
 *            on the device;
 *   Vrehihi_d has space for the highest doubles of the real parts
 *            of the Householder vectors on the device;
 *   Vrelohi_d has space for the second highest doubles of the real parts
 *            of the Householder vectors on the device;
 *   Vrehilo_d has space for the second lowest doubles of the real parts
 *            of the Householder vectors on the device;
 *   Vrelolo_d has space for the lowest doubles of the real parts
 *            of the Householder vectors on the device;
 *   Vimhihi_d has space for the highest doubles of the imaginary parts
 *            of the Householder vectors on the device;
 *   Vimlohi_d has space for the second highest doubles of the imaginary
 *            parts of the Householder vectors on the device;
 *   Vimhilo_d has space for the second lowest doubles of the imaginary parts
 *            of the Householder vectors on the device;
 *   Vimlolo_d has space for the lowest doubles of the imaginary parts
 *            of the Householder vectors on the device;
 *   betahihi_h has space for the highest doubles of the betas if verbose;
 *   betalohi_h has space for the second highest doubles of the betas,
 *            if verbose;
 *   betahilo_h has space for the second lowest doubles of the betas,
 *            if verbose;
 *   betalolo_h has space for the lowest doubles of the betas if verbose;
 *   betahihi_d has space on the device for the highest doubles of the betas;
 *   betalohi_d has space on the device for the second highest doubles
 *            of the betas;
 *   betahilo_d has space on the device for the second lowest doubles
 *            of the betas;
 *   betalolo_d has space on the device for the lowest doubles of the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vrehihi_d are the highest doubles of the real parts of the Householder
 *            vectors on the device;
 *   Vrelohi_d are the second highest doubles of the real parts of the
 *            Householder vectors on the device;
 *   Vrehilo_d are the second lowest doubles of the real parts of the
 *            Householder vectors on the device;
 *   Vrelolo_d are the lowest doubles of the real parts of the Householder
 *            vectors on the device;
 *   Vimhihi_d are the highest doubles of the imaginary parts of the
 *            Householder vectors on the device;
 *   Vimlohi_d are the second highest doubles of the imaginary parts of the
 *            Householder vectors on the device;
 *   Vimhilo_d are the second lowest doubles of the imaginary parts of the
 *            Householder vectors on the device;
 *   Vimlolo_d are the lowest doubles of the imaginary parts of the
 *            Householder vectors on the device;
 *   betahihi_h has the vector of highest doubles of the betas, if verbose;
 *   betalohi_h has the vector of second highest doubles of the betas,
 *            if verbose;
 *   betahilo_h has the vector of second lowest doubles of the betas,
 *            if verbose;
 *   betalolo_h has the vector of lowest doubles of the betas, if verbose;
 *   betahihi_d has the highest doubles of the next beta constant;
 *   betalohi_d has the second highest doubles of the next beta constant;
 *   betahilo_d has the second lowest doubles of the next beta constant;
 *   betalolo_d has the lowest doubles of the next beta constant;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl4_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *RTdotvhihi_h, double *RTdotvlohi_h,
   double *RTdotvhilo_h, double *RTdotvlolo_h,
   double *RTdotvhihi_d, double *RTdotvlohi_d,
   double *RTdotvhilo_d, double *RTdotvlolo_d,
   double *whihi_h, double *wlohi_h, double *whilo_h, double *wlolo_h,
   double *whihi_d, double *wlohi_d, double *whilo_d, double *wlolo_d,
   double *RTvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernels to update one tile, on real data.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the reduced matrix is returned on the host.
 *
 * REQUIRED : nrows - colidx > szt.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   colidx   global index of the current column;
 *   k        index of the current tile;
 *   L        local index of the column in the current tile;
 *   Ahihi_h  highest doubles of the matrix on the host;
 *   Alohi_h  second highest doubles of the matrix on the host;
 *   Ahilo_h  second lowest doubles of the matrix on the host;
 *   Alolo_h  lowest doubles of the matrix on the host;
 *   Ahihi_d  highest doubles of the matrix on the device;
 *   Alohi_d  second highest doubles of the matrix on the device;
 *   Ahilo_d  secdon lowest doubles of the matrix on the device;
 *   Alolo_d  lowest doubles of the matrix on the device;
 *   Vhihi_d  space for the highest doubles of the Householder vectors,
 *            on the device;
 *   Vlohi_d  space for the second highest doubles of the Householder
 *            vectors, on the device;
 *   Vhilo_d  space for the second lowest doubles of the Householder
 *            vectors, on the device;
 *   Vlolo_d  space for the lowest doubles of the Householder vectors,
 *            on the device;
 *   betahihi_h has space for the highest doubles of the betas if verbose;
 *   betalohi_h has space for the second highest doubles of the betas,
 *            if verbose;
 *   betahilo_h has space for the second lowest doubles of the betas,
 *            if verbose;
 *   betalolo_h has space for the lowest doubles of the betas if verbose;
 *   betahihi_d has space on the device for the highest doubles of the betas;
 *   betalohi_d has space on the device for the second highest doubles
 *            of the betas;
 *   betahilo_d has space on the device for the second lowest doubles
 *            of the betas;
 *   betalolo_d has space on the device for the lowest doubles of the betas;
 *   RTdotvhihi_h has space for the highest doubles of RTdotv,
 *            the componentwise product of R^T with v, if verbose;
 *   RTdotvlohi_h has space for the second highest doubles of RTdotv,
 *            if verbose;
 *   RTdotvhilo_h has space for the second lowest doubles of RTdotv,
 *            if verbose;
 *   RTdotvlolo_h has space for the lowest doubles of the RTdotv, if verbose;
 *   RTdotvhihi_d has space for the highest doubles of RTdotv, on the device;
 *   RTdotvlohi_d has space for the second highest doubles of RTdotv, 
 *            on the device
 *   RTdotvhilo_d has space for the second lowest doubles of RTdotv, 
 *            on the device
 *   RTdotvlolo_d has space for the lowest doubles  of RTdotv, on the device;
 *   whihi_h  space for the highest doubles of beta*R^T*v,
 *            on the host, if verbose;
 *   wlohi_h  space for the second highest doubles of beta*R^T*v,
 *            on the host, if verbose;
 *   whilo_h  space for the second lowest doubles of beta*R^T*v,
 *            on the host, if verbose;
 *   wlolo_h  space for the lowest doubles of beta*R^T*v,
 *            on the host, if verbose;
 *   whihi_d  space for the highest doubles of beta*R^T*v, plus szt padding;
 *   whilo_d  space for the second highest doubles of beta*R^T*v,
 *            plus szt padding;
 *   whilo_d  space for the second lowest doubles of beta*R^T*v,
 *            plus szt padding;
 *   wlolo_d  space for the lowest doubles of beta*R^T*v, plus szt padding;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vhi_d    contains the highest doubles of the Householder vectors V,
 *            on the device;
 *   Vlohi_d  contains the second highest doubles of V, on the device;
 *   Vhilo_d  contains the second lowest doubles of V, on the device;
 *   Vlolo_d  contains the lowest doubles of V, on the device;
 *   betahihi_h is the vector of the highest doubles of betas, if verbose;
 *   betalohi_h is the vector of the second highest doubles of betas,
 *            if verbose;
 *   betahilo_h is the vector of the second lowest doubles of betas,
 *            if verbose;
 *   betalolo_h is the vector of the lowest doubles of betas, if verbose;
 *   betahihi_d has the highest double of the next beta constant;
 *   betalohi_d has the second highest double of the next beta constant;
 *   betahilo_d has the second lowest double of the next beta constant;
 *   betalolo_d has the lowest double of the next beta constant;
 *   RTdotvhihi_h stores the highest doubles of the componentwise product
 *            of R^T with v, if verbose;
 *   RTdotvlohi_h stores the second highest doubles of the RTdotv, if verbose;
 *   RTdotvhilo_h stores the second lowest doubles of RTdotv, if verbose;
 *   RTdotvlolo_h stores the lowest doubles of RTdotv, f verbose;
 *   whihi_h  stores the highest doubles of beta*R^T*v, if verbose;
 *   wlohi_h  stores the second highest doubles of beta*R^T*v, if verbose;
 *   whilo_h  stores the second lowest doubles of beta*R^T*v, if verbose;
 *   wlolo_h  stores the lowest doubles of beta*R^T*v, if verbose;
 *   RTvlapms is the elapsed time spent to compute beta*R^T*v;
 *   redlapms is the elapsed time spent to reduce one tile;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx4_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *RHdotvrehihi_h, double *RHdotvrelohi_h,
   double *RHdotvrehilo_h, double *RHdotvrelolo_h,
   double *RHdotvimhihi_h, double *RHdotvimlohi_h,
   double *RHdotvimhilo_h, double *RHdotvimlolo_h,
   double *RHdotvrehihi_d, double *RHdotvrelohi_d,
   double *RHdotvrehilo_d, double *RHdotvrelolo_d,
   double *RHdotvimhihi_d, double *RHdotvimlohi_d,
   double *RHdotvimhilo_d, double *RHdotvimlolo_d,
   double *wrehihi_h, double *wrelohi_h, double *wrehilo_h, double *wrelolo_h,
   double *wimhihi_h, double *wimlohi_h, double *wimhilo_h, double *wimlolo_h,
   double *wrehihi_d, double *wrelohi_d, double *wrehilo_d, double *wrelolo_d,
   double *wimhihi_d, double *wimlohi_d, double *wimhilo_d, double *wimlolo_d,
   double *RHvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernels to update one tile, on complex data.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the reduced matrix is returned on the host.
 *
 * REQUIRED : nrows - colidx > szt.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   colidx   global index of the current column;
 *   k        index of the current tile;
 *   L        local index of the column in the current tile;
 *   Arehihi_h are the highest doubles of the real parts of A on the host;
 *   Arelohi_h are the second highest doubles of the real parts of A
 *            on the host;
 *   Arehilo_h are the second lowest doubles of the real parts of A
 *            on the host;
 *   Arelolo_h are the lowest doubles of the real parts of A on the host;
 *   Aimhihi_h are the highest doubles of the imaginary parts of A
 *            on the host;
 *   Aimlohi_h are the second highest doubles of the imaginary parts of A
 *            on the host;
 *   Aimhilo_h are the second lowest doubles of the imaginary parts of A
 *            on the host;
 *   Aimlolo_h are the lowest doubles of the imaginary parts of A on the host;
 *   Arehihi_d are the highest doubles of the real parts of A on the device;
 *   Arelohi_d are the second highest doubles of the real parts of A
 *            on the device;
 *   Arehilo_d are the second lowest doubles of the real parts of A
 *            on the device;
 *   Arelolo_d are the lowest doubles of the real parts of A on the device;
 *   Aimhihi_d are the highest doubles of the imaginary parts of A
 *            on the device;
 *   Aimlohi_d are the second highest doubles of the imaginary parts of A,
 *            on the device;
 *   Aimhilo_d are the second lowest doubles of the imaginary parts of A,
 *            on the device;
 *   Aimlolo_d are the lowest doubles of the imaginary parts of A
 *            on the device;
 *   Vrehihi_d has space for the highest doubles of the real parts
 *            of the Householder vectors, on the device;
 *   Vrelohi_d has space for the second highest doubles of the real parts
 *            of the Householder vectors, on the device;
 *   Vrehilo_d has space for the second lowest doubles of the real parts
 *            of the Householder vectors, on the device;
 *   Vrelolo_d has space for the lowest doubles of the real parts
 *            of the Householder vectors, on the device;
 *   Vimhihi_d has space for the highest doubles of the imaginary parts
 *            of the Householder vectors, on the device;
 *   Vimlohi_d has space for the second highest doubles of the imaginary parts
 *            of the Householder vectors, on the device;
 *   Vrehilo_d has space for the second lowest doubles of the imaginary parts
 *            of the Householder vectors, on the device;
 *   Vrelolo_d has space for the lowest doubles of the imaginary parts
 *            of the Householder vectors, on the device;
 *   betahihi_h has space for the highest doubles of the betas if verbose;
 *   betalohi_h has space for the second highest doubles of the betas
 *            if verbose;
 *   betahilo_h has space for the second lowest doubles of the betas
 *            if verbose;
 *   betalolo_h has space for the lowest doubles of the betas if verbose;
 *   betahihi_d has space on the device for the highest doubles of the betas;
 *   betalohi_d has space on the device for the second highest doubles
 *            of the betas;
 *   betahilo_d has space on the device for the second lowest doubles
 *            of the betas;
 *   betalolo_d has space on the device for the lowest doubles of the betas;
 *   RHdotvrehi_h has space for the high doubles of the real parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvrelo_h has space for the low doubles of the real parts of
 *            the componentwise product of R^H with v, if verbose;
 *   RHdotvimhi_h has space for the high doubles of the imaginary parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvimlo_h has space for the low doubles of the imaginary parts of
 *            the componentwise product of R^H with v, if verbose;
 *   RHdotvrehi_d has space for high doubles of the real parts
 *            of the componentwise product of R^H with v, on the device;
 *   RHdotvrelo_d has space for the low doubles of the real parts
 *            of the componentwise product of R^H with v, on the device;
 *   RHdotvimhi_d has space for high doubles of the imaginary parts
 *            of the componentwise product of R^H with v, on the device;
 *   RHdotvimlo_d has space for the low doubles of the imaginary parts
 *            of the componentwise product of R^H with v, on the device;
 *   wrehihi_h has space for the highest doubles of the real parts
 *            of beta*R^H*v, on the host, if verbose;
 *   wrelohi_h has space for the second highest doubles of the real parts
 *            of beta*R^H*v, on the host, if verbose;
 *   wrehilo_h has space for the second lowest doubles of the real parts
 *            of beta*R^H*v on the host, if verbose;
 *   wrelolo_h has space for the lowest doubles of the real parts
 *            of beta*R^H*v on the host, if verbose;
 *   wimhihi_h has space for the highest doubles of the imaginary parts
 *            of beta*R^H*v, on the host, if verbose;
 *   wimlohi_h has space for the second highest doubles of the imaginary parts
 *            of beta*R^H*v, on the host, if verbose;
 *   wimhilo_h has space for the second lowest doubles of the imaginary parts
 *            of beta*R^H*v on the host, if verbose;
 *   wimlolo_h has space for the lowest doubles of the imaginary parts
 *            of beta*R^H*v on the host, if verbose;
 *   wrehihi_d has space for the highest doubles of the real parts
 *            of beta*R^H*v, plus szt padding;
 *   wrelohi_d has space for the second highest doubles of the real parts
 *            of beta*R^H*v, plus szt padding;
 *   wrehilo_d has space for the second lowest doubles of the real parts
 *            of beta*R^H*v, plus szt padding;
 *   wrelolo_d has space for the lowest doubles of the real parts
 *            of beta*R^H*v, plus szt padding;
 *   wimhihi_d has space for the highest doubles of the imaginary parts
 *            of beta*R^H*v, plus szt padding;
 *   wimlohi_d has space for the second highest doubles of the imaginary parts
 *            of beta*R^H*v, plus szt padding;
 *   wimhilo_d has space for the second lowest doubles of the imaginary parts
 *            of beta*R^H*v, plus szt padding;
 *   wimlolo_d has space for the lowest doubles of the imaginary parts
 *            of beta*R^H*v, plus szt padding;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vrehihi_d are the highest doubles of the real parts
 *            of the Householder vectors, on the device;
 *   Vrelohi_d are the second highest doubles of the real parts
 *            of the Householder vectors, on the device;
 *   Vrehilo_d are the second lowest doubles of the real parts
 *            of the Householder vectors, on the device;
 *   Vrelolo_d are the lowest doubles of the real parts
 *            of the Householder vectors, on the device;
 *   Vimhihi_d are the highest doubles of the imaginary parts
 *            of the Householder vectors, on the device;
 *   Vimlohi_d are the second highest doubles of the imaginary parts
 *            of the Householder vectors, on the device;
 *   Vimhilo_d are the second lowest doubles of the imaginary parts
 *            of the Householder vectors, on the device;
 *   Vimlolo_d are the lowest doubles of the imaginary parts
 *            of the Householder vectors, on the device;
 *   betahihi_h is the vector of the highest doubles of betas, if verbose;
 *   betalohi_h is the vector of the second highest doubles of betas,
 *            if verbose;
 *   betahilo_h is the vector of the second lowest doubles of betas,
 *            if verbose;
 *   betalolo_h is the vector of the lowest doubles of betas, if verbose;
 *   betahihi_d has the highest double of the next beta constant;
 *   betalohi_d has the second highest double of the next beta constant;
 *   betahilo_d has the second lowest double of the next beta constant;
 *   betalolo_d has the lowest double of the next beta constant;
 *   RHdotvrehihi_h stores the highest doubles of the real parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvrelohi_h stores the second highest doubles of the real parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvrehilo_h stores the second lowest doubles of the real parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvrelolo_h stores the lowest doubles of the real parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvimhihi_h stores the highest doubles of the imaginary parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvimlohi_h stores the second highest doubles of the imaginary parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvimhilo_h stores the second lowest doubles of the imaginary parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvimlolo_h stores the lowest doubles of the imaginary parts
 *            of the componentwise product of R^H with v, if verbose;
 *   RHdotvhihi_d stores the highest doubles of R^H*v, on the device;
 *   RHdotvlohi_d stores the second highest doubles of R^H*v, on the device;
 *   RHdotvhilo_d stores the second lowest doubles of R^H*v, on the device;
 *   RHdotvlolo_d stores the lowest doubles of R^H*v, on the device;
 *   wrehihi_h stores the highest doubles of the real parts
 *            of beta*R^H*v, if verbose;
 *   wrelohi_h stores the second highest doubles of the real parts
 *            of beta*R^H*v, if verbose;
 *   wrehilo_h stores the second lowest doubles of the real parts
 *            of beta*R^H*v, if verbose;
 *   wrelolo_h stores the lowest doubles of the real parts
 *            of beta*R^H*v, if verbose;
 *   wimhihi_h stores the highest doubles of the imaginary parts
 *            of beta*R^H*v, if verbose;
 *   wimlohi_h stores the second highest doubles of the imaginary parts
 *            of beta*R^H*v, if verbose;
 *   wimhilo_h stores the second lowest doubles of the imaginary parts
 *            of beta*R^H*v, if verbose;
 *   wimlolo_h stores the lowest doubles of the imaginary parts
 *            of beta*R^H*v, if verbose;
 *   RHvlapms is the elapsed time spent to compute beta*R^H*v;
 *   redlapms is the elapsed time spent to reduce one tile;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl4_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vhihi_h, double *Vlohi_h, double *Vhilo_h, double *Vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *Whihi_h, double *Wlohi_h, double *Whilo_h, double *Wlolo_h,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *WYThihi_h, double *WYTlohi_h, double *WYThilo_h, double *WYTlolo_h,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the W in the WY representation,
 *   on real data.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Vhihi_h  highest doubles of the Householder vectors, if verbose;
 *   Vlohi_h  second highest doubles of the Householder vectors, if verbose;
 *   Vhilo_h  second lowest doubles of the Householder vectors, if verbose;
 *   Vlolo_h  lowest doubles of the Householder vectors, if verbose;
 *   Vhihi_d  highest doubles of the Householder vectors on the device;
 *   Vlohi_d  second highest doubles of the Householder vectors on the device;
 *   Vhilo_d  second lowest doubles of the Householder vectors on the device;
 *   Vlolo_d  lowest doubles of the Householder vectors on the device;
 *   Whihi_h  space for the highest doubles of W, if verbose;
 *   Wlohi_h  space for the second highest doubles of W, if verbose;
 *   Whilo_h  space for the second lowest doubles of W, if verbose;
 *   Wlolo_h  space for the lowest doubles of W, if verbose;
 *   Whihi_d  space for highest doubles of W on the device;
 *   Wlohi_d  space for second highest doubles of W on the device;
 *   Whilo_d  space for second lowest doubles of W on the device;
 *   Wlolo_d  space for lowest doubles of W on the device;
 *   WYThihi_h has space for the highest doubles of W*Y^T, if verbose;
 *   WYTlohi_h has space for the second highest doubles of W*Y^T, if verbose;
 *   WYThilo_h has space for the second lowest doubles of W*Y^T, if verbose;
 *   WYTlolo_h has space for the lowest doubles of W*Y^T, if verbose;
 *   WYThihi_d has space for the highest doubles of W*Y^T, on the device;
 *   WYTlohi_d has space for the second highest doubles of W*Y^T,
 *            on the device;
 *   WYThilo_d has space for the second lowest doubles of W*Y^T,
 *            on the device;
 *   WYTlolo_d has space for the lowest doubles of W*Y^T, on the device;
 *   betahihi_h has space for the betas if verbose;
 *   betalohi_h has space for the betas if verbose;
 *   betahilo_h has space for the betas if verbose;
 *   betalolo_h has space for the betas if verbose;
 *   betahihi_d has space on the device for the betas;
 *   betalohi_d has space on the device for the betas;
 *   betahilo_d has space on the device for the betas;
 *   betalolo_d has space on the device for the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vhihi_h  highest doubles of Y, if verbose;
 *   Vlohi_h  second highest doubles of Y, if verbose;
 *   Vhilo_h  second lowest doubles of Y, if verbose;
 *   Vlolo_h  lowest doubles of Y, if verbose;
 *   Vhihi_d  highest doubles of Y, on the device;
 *   Vlohi_d  second highest doubles of Y, on the device;
 *   Vhilo_d  second lowest doubles of Y, on the device;
 *   Vlolo_d  lowest doubles of Y, on the device;
 *   Whihi_d  highest doubles of W, on the device;
 *   Wlohi_d  second highest doubles of W, on the device;
 *   Whilo_d  second lowest doubles of W, on the device;
 *   Wlolo_d  lowest doubles of W, on the device;
 *   Whihi_h  highest doubles of W, if verbose;
 *   Wlohi_h  second highest doubles of W, if verbose;
 *   Whilo_h  second lowest doubles of W, if verbose;
 *   Wlolo_h  lowest doubles of W, if verbose;
 *   WYThihi_h has the highest doubles of W*Y^T, if verbose;
 *   WYTlohi_h has the second highest doubles of W*Y^T, if verbose;
 *   WYThilo_h has the second lowest doubles of W*Y^T, if verbose;
 *   WYTlolo_h has the lowest doubles of W*Y^T, if verbose;
 *   WYThihi_d has the highest doubles of W*Y^T, on the device;
 *   WYTlohi_d has the second highest doubles of W*Y^T, on the device;
 *   WYThilo_d has the second lowest doubles of W*Y^T, on the device;
 *   WYTlolo_d has the lowest doubles of W*Y^T, on the device;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx4_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vrehihi_h, double *Vrelohi_h, double *Vrehilo_h, double *Vrelolo_h,
   double *Vimhihi_h, double *Vimlohi_h, double *Vimhilo_h, double *Vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *Wrehihi_h, double *Wrelohi_h, double *Wrehilo_h, double *Wrelolo_h,
   double *Wimhihi_h, double *Wimlohi_h, double *Wimhilo_h, double *Wimlolo_h,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *WYHrehihi_h, double *WYHrelohi_h,
   double *WYHrehilo_h, double *WYHrelolo_h,
   double *WYHimhihi_h, double *WYHimlohi_h,
   double *WYHimhilo_h, double *WYHimlolo_h,
   double *WYHrehihi_d, double *WYHrelohi_d,
   double *WYHrehilo_d, double *WYHrelolo_d,
   double *WYHimhihi_d, double *WYHimlohi_d,
   double *WYHimhilo_d, double *WYHimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the W in the WY representation,
 *   on complex data.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows     number of rows in the matrices V, Y, and W;
 *   ncols     equals the size of one tile, or equivalently,
 *             is the number of elements in B,
 *             and the number of columns in V, Y, and W;
 *   szt       size of one tile and the number of threads in a block;
 *   idx       index of the current tile;
 *   Vrehihi_h are the highest doubles of the real parts
 *             of the Householder vectors, if verbose;
 *   Vrelohi_h are the second highest doubles of the real parts
 *             of the Householder vectors, if verbose;
 *   Vrehilo_h are the second lowest doubles of the real parts
 *             of the Householder vectors, if verbose;
 *   Vrelolo_h are the lowest doubles of the real parts
 *             of the Householder vectors, if verbose;
 *   Vimhihi_h are the highest doubles of the imaginary parts
 *             of the Householder vectors, if verbose;
 *   Vimlohi_h are the second highest doubles of the imaginary parts
 *             of the Householder vectors, if verbose;
 *   Vimhilo_h are the second lowest doubles of the imaginary parts
 *             of the Householder vectors, if verbose;
 *   Vimlolo_h are the lowest doubles of the imaginary parts
 *             of the Householder vectors, if verbose;
 *   Vrehihi_d are the highest doubles of the the real parts of V,
 *             on the device;
 *   Vrelohi_d are the second highest doubles of the the real parts of V,
 *             on the device;
 *   Vrehilo_d are the second lowest doubles of the real parts of V,
 *             on the device;
 *   Vrelolo_d are the lowest doubles of the real parts of V,
 *             on the device;
 *   Vimhihi_d are the highest doubles of the imaginary parts of V,
 *             on the device;
 *   Vimlohi_d are the second highest doubles of the imaginary parts of V,
 *             on the device;
 *   Vimhilo_d are the second lowest doubles of the imaginary parts of V,
 *             on the device;
 *   Vimlolo_d are the lowest doubles of the imaginary parts of V,
 *             on the device;
 *   Wrehihi_h has space for the highest doubles of the real parts of W,
 *             if verbose;
 *   Wrelohi_h has space for the second highest doubles of the real parts
 *             of W, if verbose;
 *   Wrehilo_h has space for the second lowest doubles of the real parts of W,
 *             if verbose;
 *   Wrelolo_h has space for the lowest doubles of the real parts of W,
 *             if verbose;
 *   Wimhihi_h has space for the highest doubles of the imaginary parts of W,
 *             if verbose;
 *   Wimlohi_h has space for the second highest doubles of the imaginary parts
 *             of W, if verbose;
 *   Wimhilo_h has space for the second lowest doubles of the imaginary parts
 *             of W, if verbose;
 *   Wimlolo_h has space for the lowest doubles of the imaginary parts of W,
 *             if verbose;
 *   Wrehihi_d has space for the highest doubles of the real parts of W,
 *             on the device;
 *   Wrelohi_d has space for the second highest doubles of the real parts
 *             of W, on the device;
 *   Wrehilo_d has space for the second lowest doubles of the real parts of W,
 *             on the device;
 *   Wrelolo_d has space for the lowest doubles of the real parts of W,
 *             on the device;
 *   Wimhihi_d has space for the highest doubles of the imaginary parts of W,
 *            on the device;
 *   Wimlohi_d has space for the second highest doubles of the imaginary parts
 *             of W, on the device;
 *   Wimhilo_d has space for the second lowest doubles of the imaginary parts
 *             of W, on the device;
 *   Wimlolo_d has space for the lowest doubles of the imaginary parts of W,
 *             on the device;
 *   WYHrehihi_h has space for the highest doubles of the real parts
 *             of W*Y^H, if verbose;
 *   WYHrelohi_h has space for the second highest doubles of the real parts
 *             of W*Y^H, if verbose;
 *   WYHrehilo_h has space for the second lowest doubles of the real parts
 *             of W*Y^H, if verbose;
 *   WYHrelolo_h has space for the lowest doubles of the real parts
 *             of W*Y^H, if verbose;
 *   WYHimhihi_h has space for the highest doubles of the imaginary parts
 *             of W*Y^H, if verbose;
 *   WYHimlohi_h has space for the second highest doubles of the imaginary
 *             parts of W*Y^H, if verbose;
 *   WYHimhilo_h has space for the secondlowest doubles of the imaginary parts
 *             of W*Y^H, if verbose;
 *   WYHimlolo_h has space for the lowest doubles of the imaginary parts
 *             of W*Y^H, if verbose;
 *   WYHrehihi_d has space for the highest doubles of the real parts
 *             of W*Y^H, on the device;
 *   WYHrelohi_d has space for the second highest doubles of the real parts
 *             of W*Y^H, on the device;
 *   WYHrehilo_d has space for the second lowest doubles of the real parts
 *             of W*Y^H, on the device;
 *   WYHrelolo_d has space for the lowest doubles of the real parts
 *             of W*Y^H, on the device;
 *   WYHimhihi_d has space for the highest doubles of the imaginary parts
 *             of W*Y^H, on the device;
 *   WYHimlohi_d has space for the second highest doubles of the imaginary
 *             parts of W*Y^H, on the device;
 *   WYHimhilo_d has space for the second lowest doubles of the imaginary
 *             parts of W*Y^H, on the device;
 *   WYHimlolo_d has space for the lowest doubles of the imaginary parts
 *             of W*Y^H, on the device;
 *   betahihi_h has space for the betas if verbose;
 *   betalohi_h has space for the betas if verbose;
 *   betahilo_h has space for the betas if verbose;
 *   betalolo_h has space for the betas if verbose;
 *   betahihi_d has space on the device for the betas;
 *   betalohi_d has space on the device for the betas;
 *   betahilo_d has space on the device for the betas;
 *   betalolo_d has space on the device for the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vrehi_h  high doubles of the real parts of the Y matrix, if verbose;
 *   Vrelo_h  low doubles of the real parts of the Y matrix, if verbose;
 *   Vimhi_h  high doubles of the imaginary parts of the Y matrix, if verbose;
 *   Vimlo_h  low doubles of the imaginary parts of the Y matrix, if verbose;
 *   Vrehi_d  high doubles of the real parts of the Y matrix on the device;
 *   Vrelo_d  low doubles of the real parts of the Y matrix on the device;
 *   Vimhi_d  high doubles of the imaginary parts of the Y matrix,
 *            on the device;
 *   Vimlo_d  low doubles of the imaginary parts of the Y matrix,
 *            on the device;
 *   Wrehi_d  high doubles of the the real parts of W, on the device;
 *   Wrelo_d  low doubles of the the real parts of W, on the device;
 *   Wimhi_d  high doubles of the imaginary parts of W, on the device;
 *   Wimlo_d  low doubles of the imaginary parts of W, on the device;
 *   Wrehi_h  high doubles of the real parts of W, if verbose;
 *   Wrelo_h  low doubles of the real parts of W, if verbose;
 *   Wimhi_h  high doubles of the imaginary parts of W, if verbose;
 *   Wimlo_h  low doubles of the imaginary parts of W, if verbose;
 *   WYHrehi_h has the high doubles of the real parts of W*Y^H,
 *            if verbose;
 *   WYHrelo_h has the low doubles of the real parts of W*Y^H,
 *            if verbose;
 *   WYHimhi_h has the high doubles of the imaginary parts of W*Y^H,
 *            if verbose;
 *   WYHimlo_h has the low doubles of the imaginary parts of W*Y^H,
 *            if verbose;
 *   WYHrehi_d has the high doubles of the real parts of W*Y^H,
 *            on the device;
 *   WYHrelo_d has the low doubles of the real parts of W*Y^H,
 *            on the device;
 *   WYHimhi_d has the high doubles of the imaginary parts of W*Y^H,
 *            on the device;
 *   WYHimlo_d has the low doubles of the imaginary parts of W*Y^H,
 *            on the device;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl4_small_WYT
 ( int nrows, int szt,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *Yhihi_d, double *Ylohi_d, double *Yhilo_d, double *Ylolo_d,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *WYThihi_h, double *WYTlohi_h, double *WYThilo_h, double *WYTlolo_h,
   double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^T.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^T matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   Whihi_d  highest doubles of W in the WY representation;
 *   Wlohi_d  second highest doubles of W in the WY representation;
 *   Whilo_d  second lowest doubles of W in the WY representation;
 *   Wlolo_d  lowest doubles of W in the WY representation;
 *   Yhihi_d  highest doubles of the matrix Y of Householder vectors;
 *   Ylohi_d  second highest doubles of Y;
 *   Yhilo_d  second lowest doubles of Y;
 *   Ylolo_d  lowest doubles of Y;
 *   WYThihi_d has space for an nrows-by-nrows matrix on the device;
 *   WYTlohi_d has space for an nrows-by-nrows matrix on the device;
 *   WYThilo_d has space for an nrows-by-nrows matrix on the device;
 *   WYTlolo_d has space for an nrows-by-nrows matrix on the device;
 *   WYThihi_h has space for an nrows-by-nrows matrix on the host, if verbose;
 *   WYTlohi_h has space for an nrows-by-nrows matrix on the host, if verbose;
 *   WYThilo_h has space for an nrows-by-nrows matrix on the host, if verbose;
 *   WYTlolo_h has space for an nrows-by-nrows matrix on the host, if verbose;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   WYThihi_d has the highest doubles of W*Y^T on the device;
 *   WYTlohi_d has the highest doubles of W*Y^T on the device;
 *   WYThilo_d has the lowest doubles of W*Y^T on the device;
 *   WYTlolo_d has the lowest doubles of W*Y^T on the device;
 *   WYThihi_h has the highest doubles of W*Y^T, if verbose;
 *   WYTlohi_h has the highest doubles of W*Y^T, if verbose;
 *   WYThilo_h has the lowest doubles of W*Y^T, if verbose;
 *   WYTlolo_h has the lowest doubles of W*Y^T, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_cmplx4_small_WYH
 ( int nrows, int szt,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *Yrehihi_d, double *Yrelohi_d, double *Yrehilo_d, double *Yrelolo_d,
   double *Yimhihi_d, double *Yimlohi_d, double *Yimhilo_d, double *Yimlolo_d,
   double *WYTrehihi_d, double *WYTrelohi_d,
   double *WYTrehilo_d, double *WYTrelolo_d,
   double *WYTimhihi_d, double *WYTimlohi_d,
   double *WYTimhilo_d, double *WYTimlolo_d,
   double *WYTrehihi_h, double *WYTrelohi_h,
   double *WYTrehilo_h, double *WYTrelolo_h,
   double *WYTimhihi_h, double *WYTimlohi_h,
   double *WYTimhilo_h, double *WYTimlolo_h,
   double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^H.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^H matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   Wrehihi_d has the highest doubles of the real parts of W;
 *   Wrelohi_d has the second highest doubles of the real parts of W;
 *   Wrehilo_d has the second lowest doubles of the real parts of W;
 *   Wrelolo_d has the lowest doubles of the real parts of W;
 *   Wimhihi_d has the highest doubles of the imaginary parts of W;
 *   Wimlohi_d has the second highest doubles of the imaginary parts of W;
 *   Wimhilo_d has the second lowest doubles of the imaginary parts of W;
 *   Wimlolo_d has the lowest doubles of the imaginary parts of W;
 *   Yrehi_d  has the highest doubles of the real parts of the matrix Y
 *            of the Householder vectors;
 *   Yrelohi_d has the second highest doubles of the real parts of Y;
 *   Yrehilo_d has the second lowest doubles of the real parts of Y;
 *   Yrelolo_d has the lowest doubles of the real parts of Y;
 *   Yimhihi_d has the highest doubles of the imaginary parts of  Y;
 *   Yimlohi_d has the second highest doubles of the imaginary parts of  Y;
 *   Yimhilo_d has the second lowest doubles of the imaginary parts of  Y;
 *   Yimlolo_d has the lowest doubles of the imaginary parts of  Y;
 *   WYTrehihi_d has space for an nrows-by-nrows matrix on the device;
 *   WYTrelohi_d has space for an nrows-by-nrows matrix on the device;
 *   WYTrehilo_d has space for an nrows-by-nrows matrix on the device;
 *   WYTrelolo_d has space for an nrows-by-nrows matrix on the device;
 *   WYTimhihi_d has space for an nrows-by-nrows matrix on the device;
 *   WYTimlohi_d has space for an nrows-by-nrows matrix on the device;
 *   WYTimhilo_d has space for an nrows-by-nrows matrix on the device;
 *   WYTimlolo_d has space for an nrows-by-nrows matrix on the device;
 *   WYTrehihi_h has space for an nrows-by-nrows matrix on the host,
 *            if verbose;
 *   WYTrelohi_h has space for an nrows-by-nrows matrix on the host;
 *   WYTrehilo_h has space for an nrows-by-nrows matrix on the host;
 *   WYTrelolo_h has space for an nrows-by-nrows matrix on the host;
 *   WYTimhihi_h has space for an nrows-by-nrows matrix on the host;
 *   WYTimlohi_h has space for an nrows-by-nrows matrix on the host;
 *   WYTimhilo_h has space for an nrows-by-nrows matrix on the host;
 *   WYTimlolo_h has space for an nrows-by-nrows matrix on the host;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   WYTrehihi_d are the highest doubles of the real parts of W*Y^T
 *            on the device;
 *   WYTrelohi_d are the second highest doubles of the real parts of W*Y^T
 *            on the device;
 *   WYTrehilo_d are the second lowest doubles of the real parts of W*Y^T
 *            on the device;
 *   WYTrelolo_d are the lowest doubles of the real parts of W*Y^T
 *            on the device;
 *   WYTimhihi_d are the highest doubles of the imaginary parts
 *            of W*Y^T on the device;
 *   WYTimlohi_d are the second highest doubles of the imaginary parts
 *            of W*Y^T on the device;
 *   WYTimhilo_d are the second lowest doubles of the imaginary parts
 *            of W*Y^T on the device;
 *   WYTimlolo_d are the lowest doubles of the imaginary parts
 *            of W*Y^T on the device;
 *   WYTrehihi_h are the highest doubles of the real parts of W*Y^T,
 *            if verbose;
 *   WYTrelohi_h are the second highest doubles of the real parts of W*Y^T,
 *            if verbose;
 *   WYTrehilo_h are the second lowest doubles of the real parts of W*Y^T,
 *            if verbose;
 *   WYTrelolo_h are the lowest doubles of the real parts of W*Y^T,
 *            if verbose;
 *   WYTimhihi_h are the highest doubles of the imaginary parts
 *            of W*Y^T, if verbose;
 *   WYTimlohi_h are the second highest doubles of the imaginary parts
 *            of W*Y^T, if verbose;
 *   WYTimhilo_h are the second lowest doubles of the imaginary parts
 *            of W*Y^T, if verbose;
 *   WYTimlolo_h are the lowest doubles of the imaginary parts
 *            of W*Y^T, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl4_small_YWT
 ( int nrows, int szt, int idx,
   double *Yhihi_d, double *Ylohi_d, double *Yhilo_d, double *Ylolo_d,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *YWThihi_d, double *YWTlohi_d, double *YWThilo_d, double *YWTlolo_d,
   double *YWThihi_h, double *YWTlohi_h, double *YWThilo_h, double *YWTlolo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^T.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^T matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Yhihi_d  highest doubles of the matrix Y of Householder vectors;
 *   Ylohi_d  second highest doubles of Y;
 *   Yhilo_d  second lowest doubles of Y;
 *   Ylolo_d  lowest doubles of Y;
 *   Whihi_d  highest doubles of W in the WY representation;
 *   Wlohi_d  second highest doubles of W;
 *   Whilo_d  second lowest doubles of W;
 *   Wlolo_d  lowest doubles of W;
 *   YWThihi_d has space for an nrows-by-nrows matrix on the device;
 *   YWTlohi_d has space for an nrows-by-nrows matrix on the device;
 *   YWThilo_d has space for an nrows-by-nrows matrix on the device;
 *   YWTlolo_d has space for an nrows-by-nrows matrix on the device;
 *   YWThihi_h has space for an nrows-by-nrows matrix on the host, if verbose;
 *   YWTlohi_h has space for an nrows-by-nrows matrix on the host, if verbose;
 *   YWThilo_h has space for an nrows-by-nrows matrix on the host, if verbose;
 *   YWTlolo_h has space for an nrows-by-nrows matrix on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWThihi_d are the highest doubles of Y*W^T on the device;
 *   YWTlohi_d are the second highest doubles of Y*W^T on the device;
 *   YWThilo_d are the second lowest doubles of Y*W^T on the device;
 *   YWTlolo_d are the lowest doubles of Y*W^T on the device;
 *   YWThihi_h are the highest doubles of Y*W^T, if verbose;
 *   YWTlohi_h are the second highest doubles of Y*W^T, if verbose;
 *   YWThilo_h are the second lowest doubles of Y*W^T, if verbose;
 *   YWTlolo_h are the lowest doubles of Y*W^T, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx4_small_YWH
 ( int nrows, int szt, int idx,
   double *Yrehihi_d, double *Yrelohi_d, double *Yrehilo_d, double *Yrelolo_d,
   double *Yimhihi_d, double *Yimlohi_d, double *Yimhilo_d, double *Yimlolo_d,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *YWTrehihi_d, double *YWTrelohi_d,
   double *YWTrehilo_d, double *YWTrelolo_d,
   double *YWTimhihi_d, double *YWTimlohi_d,
   double *YWTimhilo_d, double *YWTimlolo_d,
   double *YWTrehihi_h, double *YWTrelohi_h,
   double *YWTrehilo_h, double *YWTrelolo_h,
   double *YWTimhihi_h, double *YWTimlohi_h,
   double *YWTimhilo_h, double *YWTimlolo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^H.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^H matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Yrehihi_d has the highest doubles of the real parts of Y;
 *   Yrelohi_d has the second highest doubles of the real parts of Y;
 *   Yrehilo_d has the second lowest doubles of the real parts of Y;
 *   Yrelolo_d has the lowest doubles of the real parts of Y;
 *   Yimhihi_d has the highest doubles of the imaginary parts of Y;
 *   Yimlohi_d has the second highest doubles of the imaginary parts of Y;
 *   Yimhilo_d has the second lowest doubles of the imaginary parts of Y;
 *   Yimlolo_d has the lowest doubles of the imaginary parts of Y;
 *   Wrehihi_d has the highest doubles of the real parts of W, on the host;
 *   Wrelohi_d has the second highest doubles of the real parts of W; 
 *   Wrehilo_d has the second lowest doubles of the real parts of W; 
 *   Wrelolo_d has the lowest doubles of the real parts of W; 
 *   Wimhihi_d has the highest doubles of the imaginary parts of W,
 *            on the device;
 *   Wimlohi_d has the second highest doubles of the imaginary parts o  W;
 *   Wimhilo_d has the second lowest doubles of the imaginary parts of W;
 *   Wimlolo_d has the lowest doubles of the imaginary parts of W;
 *   YWHrehihi_d has space for an nrows-by-nrows matrix on the device;
 *   YWHrelohi_d has space for an nrows-by-nrows matrix on the device;
 *   YWHrehilo_d has space for an nrows-by-nrows matrix on the device;
 *   YWHrelolo_d has space for an nrows-by-nrows matrix on the device;
 *   YWHimhihi_d has space for an nrows-by-nrows matrix on the device;
 *   YWHimlohi_d has space for an nrows-by-nrows matrix on the device;
 *   YWHimhilo_d has space for an nrows-by-nrows matrix on the device;
 *   YWHimlolo_d has space for an nrows-by-nrows matrix on the device;
 *   YWHrehihi_h has space for an nrows-by-nrows matrix on the host,
 *            if verbose;
 *   YWHrelohi_h has space for an nrows-by-nrows matrix on the host;
 *   YWHrehilo_h has space for an nrows-by-nrows matrix on the host;
 *   YWHrelolo_h has space for an nrows-by-nrows matrix on the host;
 *   YWHimhihi_h has space for an nrows-by-nrows matrix on the host;
 *   YWHimlohi_h has space for an nrows-by-nrows matrix on the host;
 *   YWHimhilo_h has space for an nrows-by-nrows matrix on the host;
 *   YWHimlolo_h has space for an nrows-by-nrows matrix on the host;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWHrehihi_d are the highest doubles of the real parts of Y*W^H,
 *            on the device;
 *   YWHrelohi_d are the second highest doubles of the real parts of Y*W^H,
 *            on the device;
 *   YWHrehilo_d are the second lowest doubles of the real parts of Y*W^H,
 *            on the device;
 *   YWHrelolo_d are the lowest doubles of the real parts of Y*W^H,
 *            on the device;
 *   YWHimhihi_d are the highest doubles of imaginary parts of Y*W^H,
 *            on the device;
 *   YWHimlohi_d are the highest doubles of imaginary parts of Y*W^H,
              on the device;
 *   YWHimhilo_d are the second lowest doubles of imaginary parts of Y*W^H,
 *            on the device;
 *   YWHimlolo_d are the lowest doubles of imaginary parts of Y*W^H,
 *            on the device;
 *   YWHrehihi_h are the highest doubles of the real parts of Y*W^H,
 *            if verbose;
 *   YWHrelohi_h are the second highest doubles of the real parts of Y*W^H,
 *            if verbose;
 *   YWHrehilo_h are the second lowest doubles of the real parts of Y*W^H,
 *            if verbose;
 *   YWHrelolo_h are the lowest doubles of the real parts of Y*W^H,
 *            if verbose;
 *   YWHimhihi_h are the highest doubles of imaginary parts of Y*W^H,
 *            if verbose;
 *   YWHimlohi_h are the second highest doubles of imaginary parts of Y*W^H,
 *            if verbose;
 *   YWHimhilo_h are the second lowest doubles of imaginary parts of Y*W^H,
 *            if verbose;
 *   YWHimlolo_h are the lowest doubles of imaginary parts of Y*W^H,
 *            if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl4_small_QWYT
 ( int dim, int szt, int idx,
   double *Qhihi_d, double *Qlohi_d, double *Qhilo_d, double *Qlolo_d,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *QWYThihi_d, double *QWYTlohi_d,
   double *QWYThilo_d, double *QWYTlolo_d,
   double *QWYThihi_h, double *QWYTlohi_h,
   double *QWYThilo_h, double *QWYTlolo_h,
   double *Qhihi_h, double *Qlohi_h, double *Qhilo_h, double *Qlolo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute Q*WYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the Q*WYT matrix is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q and WYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Qhihi_d  highest doubles of Q, a dim-by-dim matrix, on the device;
 *   Qlohi_d  second highest doubles of Q, on the device;
 *   Qhilo_d  second lowest doubles of Q, on the device;
 *   Qlolo_d  lowest doubles of Q, on the device;
 *   WYThihi_d are the highest doubles of W*Y^T, on the device;
 *   WYTlohi_d are the second highest doubles of W*Y^T, on the device;
 *   WYThilo_d are the second lowest doubles of W*Y^T, on the device;
 *   WYTlolo_d are the lowest doubles of W*Y^T, on the device;
 *   QWYThi_d has space for Q*WYT, on the device;
 *   QWYTlo_d has space for Q*WYT, on the device;
 *   QWYThi_h has space for Q*WYT, on the host, if verbose;
 *   QWYTlo_h has space for Q*WYT, on the host, if verbose;
 *   Qhihi_h  if verbose, then used to print Q before the product;
 *   Qlohi_h  if verbose, then used to print Q before the product;
 *   Qhilo_h  if verbose, then used to print Q before the product;
 *   Qlolo_h  if verbose, then used to print Q before the product;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   QWYThihi_d has the highest doubles of Q*WYT on the device;
 *   QWYTlohi_d has the second highest doubles of Q*WYT on the device;
 *   QWYThilo_d has the second lowest doubles of Q*WYT on the device;
 *   QWYTlolo_d has the lowest doubles of Q*WYT on the device;
 *   QWYThihi_h has the highest doubles of Q*WYT, if verbose;
 *   QWYTlohi_h has the second highest doubles of Q*WYT, if verbose;
 *   QWYThilo_h has the second lowest doubles of Q*WYT, if verbose;
 *   QWYTlolo_h has the lowest doubles of Q*WYT, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx4_small_QWYH
 ( int dim, int szt, int idx,
   double *Qrehihi_d, double *Qrelohi_d, double *Qrehilo_d, double *Qrelolo_d,
   double *Qimhihi_d, double *Qimlohi_d, double *Qimhilo_d, double *Qimlolo_d,
   double *WYTrehihi_d, double *WYTrelohi_d,
   double *WYTrehilo_d, double *WYTrelolo_d,
   double *WYTimhihi_d, double *WYTimlohi_d,
   double *WYTimhilo_d, double *WYTimlolo_d,
   double *QWYTrehihi_d, double *QWYTrelohi_d,
   double *QWYTrehilo_d, double *QWYTrelolo_d,
   double *QWYTimhihi_d, double *QWYTimlohi_d,
   double *QWYTimhilo_d, double *QWYTimlolo_d,
   double *QWYTrehihi_h, double *QWYTrelohi_h,
   double *QWYTrehilo_h, double *QWYTrelolo_h,
   double *QWYTimhihi_h, double *QWYTimlohi_h,
   double *QWYTimhilo_h, double *QWYTimlolo_h,
   double *Qrehihi_h, double *Qrelohi_h, double *Qrehilo_h, double *Qrelolo_h,
   double *Qimhihi_h, double *Qimlohi_h, double *Qimhilo_h, double *Qimlolo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute Q*WYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the Q*WYT matrix is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q and WYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Qrehihi_d are the highest doubles of the real parts of Q, on the device;
 *   Qrelohi_d are the second highest doubles of the real parts of Q;
 *   Qrehilo_d are the second lowest doubles of the real parts of Q;
 *   Qrelolo_d are the lowest doubles of the real parts of Q;
 *   Qimhihi_d are the highest doubles of the imaginary parts of Q;
 *   Qimlohi_d are the second highest doubles of the imaginary parts of Q;
 *   Qimhilo_d are the second lowest doubles of the imaginary parts of Q;
 *   Qimlolo_d are the lowest doubles of the imaginary parts of Q;
 *   WYTrehihi_d are the highest doubles of the real parts of W*Y^T,
 *            on the device;
 *   WYTrelohi_d are the second highest doubles of real parts of W*Y^T;
 *   WYTrehilo_d are the second lowest doubles of real parts of W*Y^T;
 *   WYTrelolo_d are the lowest doubles of real parts of W*Y^T;
 *   WYTimhihi_d are the highest doubles of the imaginary parts of W*Y^T;
 *   WYTimlohi_d are the second highest doubles of the imaginary parts
 *            of W*Y^T;
 *   WYTimhilo_d are the second lowest doubles of the imaginary parts
 *            of W*Y^T;
 *   WYTimlolo_d are the lowest doubles of the imaginary parts of W*Y^T;
 *   QWYTrehi_d has space for the highest doubles of the real
 *            parts for Q*WYT, on the device;
 *   QWYTrehi_d has space for the second highest doubles of the real
 *            parts for Q*WYT, on the device;
 *   QWYTrelo_d has space for the second lowest doubles of the real
 *            parts for Q*WYT, on the device;
 *   QWYTrelo_d has space for the lowest doubles of the real
 *            parts for Q*WYT, on the device;
 *   QWYTimhihi_d has space for the highest doubles of the imaginary
 *            parts of Q*WYT, on the device;
 *   QWYTimlohi_d has space for the second highest doubles of the imaginary
 *            parts of Q*WYT, on the device;
 *   QWYTimhilo_d has space for the second lowest doubles of the imaginary
 *            parts of Q*WYT, on the device;
 *   QWYTimlolo_d has space for the lowest doubles of the imaginary
 *            parts of Q*WYT, on the device;
 *   QWYTrehihi_h has space for the highest doubles of the real
 *            parts of Q*WYT, on the host, if verbose;
 *   QWYTrelohi_h has space for the second highest doubles of the real
 *            parts of Q*WYT, on the host, if verbose;
 *   QWYTrehilo_h has space for the second lowest doubles of the real
 *            parts of Q*WYT, on the host, if verbose;
 *   QWYTrelolo_h has space for the lowest doubles of the real
 *            parts of Q*WYT, on the host, if verbose;
 *   QWYTimhihi_h has space for the highest doubles of the imaginary
 *            parts of Q*WYT, on the host, if verbose;
 *   QWYTimlohi_h has space for the second highest doubles of the imaginary
 *            parts of Q*WYT, on the host, if verbose;
 *   QWYTimhilo_h has space for the second lowest doubles of the imaginary
 *            parts of Q*WYT, on the host, if verbose;
 *   QWYTimlolo_h has space for the lowest doubles of the imaginary
 *            parts of Q*WYT, on the host, if verbose;
 *   Qrehihi_h is used to print the highest doubles
 *            of the real parts of Q, if verbose;
 *   Qrelohi_h is used to print the second highest doubles
 *            of the real parts of Q, if verbose;
 *   Qrehilo_h is used to print the second lowest doubles
 *            of the real parts of Q, if verbose;
 *   Qrelolo_h is used to print the lowest doubles
 *            of the real parts of Q, if verbose;
 *   Qimhihi_h is used to print the highest doubles
 *            of the imaginary parts of Q, if verbose;
 *   Qimlohi_h is used to print the second highest doubles
 *            of the imaginary parts of Q, if verbose;
 *   Qimhilo_h is used to print the second lowest doubles
 *            of the imaginary parts of Q, if verbose;
 *   Qimlolo_h is used to print the lowest doubles
 *            of the imaginary parts of Q, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   QWYTrehi_d are the high doubles of the real parts of Q*WYT on the device;
 *   QWYTrelo_d are the low doubles of the real parts of Q*WYT on the device;
 *   QWYTimhi_d are the high doubles of the imaginary parts of Q*WYT,
 *            on the device;
 *   QWYTimlo_d are the low doubles of the imaginary parts of Q*WYT,
 *            on the device;
 *   QWYTrehi_h are the high doubles of the real parts of Q*WYT, if verbose;
 *   QWYTrelo_h are the low doubles of the real parts of Q*WYT, if verbose;
 *   QWYTimhi_h are the high doubles of the imaginary parts of Q*WYT,
 *            if verbose;
 *   QWYTimlo_h are the low doubles of the imaginary parts of Q*WYT,
 *            if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl4_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWThihi_d, double *YWTlohi_d, double *YWThilo_d, double *YWTlolo_d,
   double *Chihi_d, double *Clohi_d, double *Chilo_d, double *Clolo_d,
   double *YWTChihi_d, double *YWTClohi_d,
   double *YWTChilo_d, double *YWTClolo_d,
   double *YWTChihi_h, double *YWTClohi_h,
   double *YWTChilo_h, double *YWTClolo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute YWT*C.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the YWT*C matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows of the matrices C, YWT, YWTC,
 *            and the number of columns of the matrix YWT;
 *   ncols    number of columns of the matrix C;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   YWThihi_d are the highest doubles of Y*W^T, on the device;
 *   YWTlohi_d are the second highest doubles of Y*W^T, on the device;
 *   YWThilo_d are the second lowest doubles of Y*W^T, on the device;
 *   YWTlolo_d are the lowest doubles of Y*W^T, on the device;
 *   Chihi_d are the highest doubles of C, an nrows-by-ncols matrix,
 *            on the device;
 *   Clohi_d are the second highest doubles of C, on the device;
 *   Chilo_d are the second lowest doubles of C, on the device;
 *   Clolo_d are hte lowest doubles of C, on the device;
 *   YWTChihi_d has space for YWT*C, on the device;
 *   YWTClohi_d has space for YWT*C, on the device;
 *   YWTChilo_d has space for YWT*C, on the device;
 *   YWTClolo_d has space for YWT*C, on the device;
 *   YWTChihi_h has space for YWT*C, on the host, if verbose;
 *   YWTClohi_h has space for YWT*C, on the host, if verbose;
 *   YWTChilo_h has space for YWT*C, on the host, if verbose;
 *   YWTClolo_h has space for YWT*C, on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWTChihi_d has the highest doubles of YWT*C on the device;
 *   YWTClohi_d has the second highest doubles of YWT*C on the device;
 *   YWTChilo_d has the second lowest doubles of YWT*C on the device;
 *   YWTClolo_d has the lowest doubles of YWT*C on the device;
 *   YWTChihi_h has the highest doubles of YWT*C, if verbose;
 *   YWTClohi_h has the second highest doubles of YWT*C, if verbose;
 *   YWTChilo_h has the second lowest doubles of YWT*C, if verbose;
 *   YWTClolo_h has the lowest doubles of YWT*C, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx4_small_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTrehihi_d, double *YWTrelohi_d,
   double *YWTrehilo_d, double *YWTrelolo_d,
   double *YWTimhihi_d, double *YWTimlohi_d,
   double *YWTimhilo_d, double *YWTimlolo_d,
   double *Crehihi_d, double *Crelohi_d, double *Crehilo_d, double *Crelolo_d,
   double *Cimhihi_d, double *Cimlohi_d, double *Cimhilo_d, double *Cimlolo_d,
   double *YWTCrehihi_d, double *YWTCrelohi_d,
   double *YWTCrehilo_d, double *YWTCrelolo_d,
   double *YWTCimhihi_d, double *YWTCimlohi_d,
   double *YWTCimhilo_d, double *YWTCimlolo_d,
   double *YWTCrehihi_h, double *YWTCrelohi_h,
   double *YWTCrehilo_h, double *YWTCrelolo_h,
   double *YWTCimhihi_h, double *YWTCimlohi_h,
   double *YWTCimhilo_h, double *YWTCimlolo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute YWT*C.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the YWT*C matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows of the matrices C, YWT, YWTC,
 *            and the number of columns of the matrix YWT;
 *   ncols    number of columns of the matrix C;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   YWTrehi_d are the high doubles of the real parts of Y*W^T, on the device;
 *   YWTrelo_d are the low doubles of the real parts of Y*W^T, on the device;
 *   YWTimhi_d are the high doubles of the imaginary parts of Y*W^T,
 *            on the device;
 *   YWTimlo_d are the low doubles of the imaginary parts of Y*W^T,
 *            on the device;
 *   Crehi_d  high doubles of the real parts of an nrows-by-ncols matrix,
 *            on the device;
 *   Crelo_d  low doubles of the real parts of an nrows-by-ncols matrix,
 *            on the device;
 *   Cimhi_d  high doubles of the imaginary parts of an nrows-by-ncols matrix,
 *            on the device;
 *   Cimlo_d  low doubles of the imaginary parts of an nrows-by-ncols matrix,
 *            on the device;
 *   YWTCrehi_d has space for the high doubles of the real parts of YWT*C,
 *            on the device;
 *   YWTCrelo_d has space for the low doubles of the real parts of YWT*C,
 *            on the device;
 *   YWTCimhi_d has space for the high doubles of the imaginary parts of YWT*C,
 *            on the device;
 *   YWTCimlo_d has space for the low doubles of the imaginary parts of YWT*C,
 *            on the device;
 *   YWTCrehi_h has space for the high doubles of the real parts of YWT*C,
 *            on the host, if verbose;
 *   YWTCrelo_h has space for the low doubles of the real parts of YWT*C,
 *            on the host, if verbose;
 *   YWTCimhi_h has space for the high doubles of the imaginary parts of YWT*C,
 *            on the host, if verbose;
 *   YWTCimlo_h has space for the low doubles of the imaginary parts of YWT*C,
 *            on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWTCrehihi_d are the highest doubles of the real parts of YWT*C
 *            on the device;
 *   YWTCrelohi_d are the second highest doubles of the real parts of YWT*C;
 *   YWTCrehilo_d are the second lowest doubles of the real parts of YWT*C;
 *   YWTCrelolo_d are the lowest doubles of the real parts of YWT*C;
 *   YWTCimhihi_d are the highest doubles of the imaginary parts of YWT*C;
 *   YWTCimlohi_d are the second highest doubles of the imaginary parts
 *            of YWT*C;
 *   YWTCimhilo_d are the second lowest doubles of the imaginary parts
 *            of YWT*C;
 *   YWTCimlolo_d are the lowest doubles of the imaginary parts of YWT*C;
 *   YWTCrehihi_h are the highest doubles of the real parts of YWT*C,
 *            if verbose;
 *   YWTCrelohi_h are the second highest doubles of the real parts of YWT*C,
 *            if verbose;
 *   YWTCrehilo_h are the second lowest doubles of the real parts of YWT*C,
 *            if verbose;
 *   YWTCrelolo_h are the lowest doubles of the real parts of YWT*C,
 *            if verbose;
 *   YWTCimhihi_h are the highest doubles of the imaginary parts of YWT*C,
 *            if verbose;
 *   YWTCimlohi_h are the second highest doubles of the imaginary parts
 *            of YWT*C, if verbose;
 *   YWTCimhilo_h are the second lowest doubles of the imaginary parts
 *            of YWT*C, if verbose;
 *   YWTCimlo_h are the lowest doubles of the imaginary parts of YWT*C,
 *            if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl4_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qhihi_d, double *Qlohi_d, double *Qhilo_d, double *Qlolo_d,
   double *QWYThihi_d, double *QWYTlohi_d,
   double *QWYThilo_d, double *QWYTlolo_d,
   double *Qhihi_h, double *Qlohi_h, double *Qhilo_h, double *Qlolo_h,
   double *lapms, long long int *add, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update Q as Q + QWYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated Q is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q,
 *            number of rows in QWYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Qhihi_d  dim-by-dim matrix, on the device;
 *   Qlohi_d  dim-by-dim matrix, on the device;
 *   Qhilo_d  dim-by-dim matrix, on the device;
 *   Qlolo_d  dim-by-dim matrix, on the device;
 *   QWYThihi_d has the highest doubles of Q*W*Y^T, on the device;
 *   QWYTlohi_d has the second highest doubles of Q*W*Y^T, on the device;
 *   QWYThilo_d has the second lowest doubles of Q*W*Y^T, on the device;
 *   QWYTlolo_d has the lowest doubles of Q*W*Y^T, on the device;
 *   add      current number of additions and subtractions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhihi_d  highest doubles of the updated Q on the device;
 *   Qlohi_d  second highest doubles of the updated Q on the device;
 *   Qhilo_d  second lowest doubles of the updated Q on the device;
 *   Qlolo_d  lowest doubles of the updated Q on the device;
 *   Qhihi_h  highest doubles of the updated Q, if verbose;
 *   Qlohi_h  second highest doubles of the updated Q, if verbose;
 *   Qhilo_h  second lowest doubles of the updated Q, if verbose;
 *   Qlolo_h  lowest doubles of the updated Q, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions. */

void GPU_cmplx4_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qrehihi_d, double *Qrelohi_d, double *Qrehilo_d, double *Qrelolo_d,
   double *Qimhihi_d, double *Qimlohi_d, double *Qimhilo_d, double *Qimlolo_d,
   double *QWYTrehihi_d, double *QWYTrelohi_d,
   double *QWYTrehilo_d, double *QWYTrelolo_d,
   double *QWYTimhihi_d, double *QWYTimlohi_d,
   double *QWYTimhilo_d, double *QWYTimlolo_d,
   double *Qrehihi_h, double *Qrelohi_h, double *Qrehilo_h, double *Qrelolo_h,
   double *Qimhihi_h, double *Qimlohi_h, double *Qimhilo_h, double *Qimlolo_h,
   double *lapms, long long int *add, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update Q as Q + QWYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated Q is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q,
 *            number of rows in QWYT;
 *   rowdim   number of columns in QWYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Qrehihi_d has the highest doubles of the real parts of Q, on the device;
 *   Qrelohi_d has the second highest doubles of the real parts of Q;
 *   Qrehilo_d has the second lowest doubles of the real parts of Q;
 *   Qrelolo_d has the lowest doubles of the real parts of Q, on the device;
 *   Qimhihi_d has the highest doubles of the imaginary parts of Q;
 *   Qimlohi_d has the second highest doubles of the imaginary parts of Q;
 *   Qimhilo_d has the second lowest doubles of the imaginary parts of Q;
 *   Qimlolo_d has the lowest doubles of the imaginary parts of Q;
 *   QWYTrehihi_d are the highest doubles of the real parts of Q*W*Y^T,
 *            on the device;
 *   QWYTrelohi_d are the second highest doubles of the real parts of Q*W*Y^T,
 *            on the device;
 *   QWYTrehilo_d are the second lowest doubles of the real parts of Q*W*Y^T,
 *            on the device;
 *   QWYTrelolo_d are the lowest doubles of the real parts of Q*W*Y^T,
 *            on the device;
 *   QWYTimhihi_d are the highest doubles of the imaginary parts of Q*W*Y^T,
 *            on the device;
 *   QWYTimlohi_d are the second highest doubles of the imaginary parts
 *            of Q*W*Y^T, on the device;
 *   QWYTimhilo_d are the second lowest doubles of the imaginary parts
 *            of Q*W*Y^T, on the device;
 *   QWYTimlolo_d are the lowest doubles of the imaginary parts of Q*W*Y^T,
 *            on the device;
 *   add      current number of additions and subtractions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qrehihi_d are the highest doubles of the real parts
 *            of the updated Q on the device;
 *   Qrelohi_d are the second highest doubles of the real parts
 *            of the updated Q on the device;
 *   Qrehilo_d are the second lowest doubles of the real parts
 *            of the updated Q on the device;
 *   Qrelolo_d are the lowest doubles of the real parts
 *            of the updated Q on the device;
 *   Qimhihi_d are the highest doubles of the imaginary parts
 *            of the updated Q on the device;
 *   Qimlohi_d are the second highest doubles of the imaginary parts
 *            of the updated Q on the device;
 *   Qimhilo_d are the second lowest doubles of the imaginary parts
 *            of the updated Q on the device;
 *   Qimlolo_d are the lowest doubles of the imaginary parts
 *            of the updated Q on the device;
 *   Qrehihi_h are the highest doubles of the real parts
 *            of the updated Q, on the host if verbose;
 *   Qrelohi_h are the second highest doubles of the real parts
 *            of the updated Q, on the host if verbose;
 *   Qrehilo_h are the second lowest doubles of the real parts
 *            of the updated Q, if verbose;
 *   Qrelolo_h are the lowest doubles of the real parts
 *            of the updated Q, if verbose;
 *   Qimhihi_h are the highest doubles of the imaginary parts
 *            of the updated Q, if verbose;
 *   Qimlohi_h are the second highest doubles of the imaginary parts
 *            of the updated Q, if verbose;
 *   Qimhilo_h are the second lowest doubles of the imaginary parts
 *            of the updated Q, if verbose;
 *   Qimlolo_h are the lowest doubles of the imaginary parts
 *            of the updated Q, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions. */

void GPU_dbl4_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rhihi_d, double *Rlohi_d, double *Rhilo_d, double *Rlolo_d,
   double *YWTChihi_d, double *YWTClohi_d,
   double *YWTChilo_d, double *YWTClolo_d,
   double *Rhihi_h, double *Rlohi_h, double *Rhilo_h, double *Rlolo_h,
   double *lapms, long long int *add, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update R as R + YWTC.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated R is returned.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   ncols    total number of columns in R and YWTC;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Rhihi_d  an nrows-by-ncols matrix, on the device;
 *   Rlohi_d  an nrows-by-ncols matrix, on the device;
 *   Rhilo_d  an nrows-by-ncols matrix, on the device;
 *   Rlolo_d  an nrows-by-ncols matrix, on the device;
 *   YWTChihi_d has the highest doubles of the product Y*W^T*C,
 *            on the device;
 *   YWTClohi_d has the second highest doubles of the product Y*W^T*C,
 *            on the device;
 *   YWTChilo_d has the second lowest doubles of the product Y*W^T*C,
 *            on the device;
 *   YWTClolo_d has the lowest doubles of the product Y*W^T*C,
 *            on the device;
 *   add      current number of additions and subtractions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Rhihi_d  highest doubles of the updated R on the device;
 *   Rlohi_d  second highest doubles of the updated R on the device;
 *   Rhilo_d  second lowest doubles of the updated R on the device;
 *   Rlolo_d  lowest doubles of the updated R on the device;
 *   Rhihi_h  highest doubles of the updated R, if verbose;
 *   Rlohi_h  second highest doubles of the updated R, if verbose;
 *   Rhilo_h  second lowest doubles of the updated R, if verbose;
 *   Rlolo_h  lowest doubles of the updated R, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions. */

void GPU_cmplx4_small_R_add_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *Rrehihi_d, double *Rrelohi_d, double *Rrehilo_d, double *Rrelolo_d,
   double *Rimhihi_d, double *Rimlohi_d, double *Rimhilo_d, double *Rimlolo_d,
   double *YWTCrehihi_d, double *YWTCrelohi_d,
   double *YWTCrehilo_d, double *YWTCrelolo_d,
   double *YWTCimhihi_d, double *YWTCimlohi_d,
   double *YWTCimhilo_d, double *YWTCimlolo_d,
   double *Rrehihi_h, double *Rrelohi_h, double *Rrehilo_h, double *Rrelolo_h,
   double *Rimhihi_h, double *Rimlohi_h, double *Rimhilo_h, double *Rimlolo_h,
   double *lapms, long long int *add, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update R as R + YWTC.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated R is returned.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   ncols    total number of columns in R and YWTC;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Rrehi_d  high doubles of the real parts of R, on the device;
 *   Rrelo_d  low doubles of the real parts of R, on the device;
 *   Rimhi_d  high doubles of the imaginary parts of R, on the device;
 *   Rimlo_d  low doubles of the imaginary parts of R, on the device;
 *   YWTCrehi_d are the high doubles of the real parts of Y*W^T*C,
 *            on the device;
 *   YWTCrelo_d are the low doubles of the real parts of Y*W^T*C,
 *            on the device;
 *   YWTCimhi_d are the high doubles of the imaginary parts of Y*W^T*C,
 *            on the device;
 *   YWTCimlo_d are the low doubles of the imaginary parts of Y*W^T*C,
 *            on the device;
 *   add      current number of additions and subtractions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Rrehi_d  high doubles of the real parts of the updated R on the device;
 *   Rrelo_d  low doubles of the real parts of the updated R on the device;
 *   Rimhi_d  high doubles of the imaginary parts of the updated R,
 *            on the device;
 *   Rimlo_d  low doubles of the imaginary parts of the updated R,
 *            on the device;
 *   Rrehi_h  high doubles of the real parts of the updated R, if verbose;
 *   Rrelo_h  low doubles of the real parts of the updated R, if verbose;
 *   Rimhi_h  high doubles of the imaginary parts of the updated R,
 *            if verbose;
 *   Rimlo_h  low doubles of the imaginary parts of the updated R,
 *            if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions. */

void GPU_dbl4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *houselapms, double *RTvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies Householder transformations in a blocked manner
 *   to compute a QR decomposition of A, on real data.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   szt      size of each block;
 *   nbt      number of tiles, ncols = szt*nbt;
 *   Ahihi    highest doubles of an nrows-by-ncols matrix A,
 *            stored as nrows arrays of ncols numbers;
 *   Alohi    second highest doubles of A;
 *   Ahilo    second lowest doubles of A;
 *   Alolo    lowest doubles of A;
 *   Qhihi    space for an nrows-by-nrows matrix;
 *   Qlohi    space for an nrows-by-nrows matrix;
 *   Qhilo    space for an nrows-by-nrows matrix;
 *   Qlolo    space for an nrows-by-nrows matrix;
 *   Rhihi    space for an nrows-by-ncols matrix;
 *   Rlohi    space for an nrows-by-ncols matrix;
 *   Rhilo    space for an nrows-by-ncols matrix;
 *   Rlolo    space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Qlohi    second highest doubles of Q;
 *   Qhilo    second lowest doubles of Q;
 *   Qlolo    lowest doubles of Q;
 *   Rhihi    highest doubles of the reduced upper triangular form R of A;
 *   Rlohi    second highest doubles of R;
 *   Rhilo    second lowest doubles of R;
 *   Rlolo    lowest doubles of R;
 *   houselapms is the elapsed time spent by the kernel
 *            to compute the Householder vector and the beta;
 *   RTvlapms is the elapsed time spent to compute beta*R^T*v;
 *   tileRlapms is the elapsed time spent by the kernel
 *            to reduce one tile;
 *   vb2Wlapms is the elapsed time spent by the kernel
 *            to compute the W representation;
 *   WYTlapms is the elapsed time spent by the kernel
 *            to compute the W*Y^T product;
 *   QWYTlapms is the elapsed time spent by the kernel
 *            to compute the Q*WYT product;
 *   Qaddlapms is the elapsed time spent by the kernel
 *            to compute Q by adding the Q*W*Y^T matrix;
 *   YWTlapms is the elapsed time spent by the kernel
 *            to compute the Y*W^T product;
 *   YWTClapms is the elapsed time spent by the kernel
 *            to compute the YWT*C product;
 *   Raddlapms is the elapsed time spent by the kernel
 *            to compute R by adding the Y*W^T*C matrix;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions;
 *   sqrtcnt  counts the number of calls to sqrt(). */

void GPU_cmplx4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *houselapms, double *RHvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYHlapms, double *QWYHlapms, double *Qaddlapms,
   double *YWHlapms, double *YWHClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies Householder transformations in a blocked manner
 *   to compute a QR decomposition of A, on complex data.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   szt      size of each block;
 *   nbt      number of tiles, ncols = szt*nbt;
 *   Arehihi  highest doubles of the real parts of A,
 *            stored as nrows arrays of ncols numbers;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the real parts of A;
 *   Arelolo  lowest doubles of the real parts of A;
 *   Aimhihi  highest doubles of the imaginary parts of A;
 *   Aimlohi  second highest doubles of the imaginary parts of A;
 *   Aimhilo  second lowest doubles of the imaginary parts of A;
 *   Aimlolo  lowest doubles of the imaginary parts of A;
 *   Qrehihi  space for an nrows-by-nrows matrix;
 *   Qrelohi  space for an nrows-by-nrows matrix;
 *   Qrehilo  space for an nrows-by-nrows matrix;
 *   Qrelolo  space for an nrows-by-nrows matrix;
 *   Qimhihi  space for an nrows-by-nrows matrix;
 *   Qimlohi  space for an nrows-by-nrows matrix;
 *   Qimhilo  space for an nrows-by-nrows matrix;
 *   Qimlolo  space for an nrows-by-nrows matrix;
 *   Rrehihi  space for an nrows-by-ncols matrix;
 *   Rrelohi  space for an nrows-by-ncols matrix;
 *   Rrehilo  space for an nrows-by-ncols matrix;
 *   Rrelolo  space for an nrows-by-ncols matrix;
 *   Rimhihi  space for an nrows-by-ncols matrix;
 *   Rimlohi  space for an nrows-by-ncols matrix;
 *   Rimhilo  space for an nrows-by-ncols matrix;
 *   Rimlolo  space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of an orthogonal matrix Q,
 *            transpose(Q)*A = R;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   Rrehihi  highest doubles of the real parts of the reduced
 *            upper triangular form R of A;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R;
 *   houselapms is the elapsed time spent by the kernel
 *            to compute the Householder vector and the beta;
 *   RHvlapms is the elapsed time spent to compute beta*R^H*v;
 *   tileRlapms is the elapsed time spent by the kernel
 *            to reduce one tile;
 *   vb2Wlapms is the elapsed time spent by the kernel
 *            to compute the W representation;
 *   WYHlapms is the elapsed time spent by the kernel
 *            to compute the W*Y^H product;
 *   QWYHlapms is the elapsed time spent by the kernel
 *            to compute the Q*WYH product;
 *   Qaddlapms is the elapsed time spent by the kernel
 *            to compute Q by adding the Q*W*Y^H matrix;
 *   YWHlapms is the elapsed time spent by the kernel
 *            to compute the Y*W^H product;
 *   YWHClapms is the elapsed time spent by the kernel
 *            to compute the YWH*C product;
 *   Raddlapms is the elapsed time spent by the kernel
 *            to compute R by adding the Y*W^H*C matrix;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions;
 *   sqrtcnt  counts the number of calls to sqrt(). */

#endif
