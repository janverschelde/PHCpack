/* The file dbl8_baqr_kernels.h specifies functions for the
 * blocked accelerated QR in octo double precision. */

#ifndef __dbl8_baqr_kernels_h__
#define __dbl8_baqr_kernels_h__

/* The constants od_shmemsize and cod_shmemsize,
 * respectively for real and complex data, determine the upper bounds
 * on the size of the largest vectors the small kernels can handle.
 * The bounds were set experimentally, based on the available amount
 * of shared memory. */

#define od_shmemsize 256
#define cod_shmemsize 128

/* The constants inner_od_shmemsize and outer_od_shmemsize,
 * respectively for the many blocks of threads and the accumulating kernel,
 * determine how munch the shared memory is consumed. */

#define inner_od_shmemsize 128
#define outer_od_shmemsize 128

/* Both constants correspond to the number of threads in a block,
 * which typically are multiples of 32. */

__global__ void dbl8_small_house
 ( double *x0hihihi, double *x0lohihi, double *x0hilohi, double *x0lolohi,
   double *x0hihilo, double *x0lohilo, double *x0hilolo, double *x0lololo,
   double *x1hihihi, double *x1lohihi, double *x1hilohi, double *x1lolohi,
   double *x1hihilo, double *x1lohilo, double *x1hilolo, double *x1lololo,
   int dim, int dimLog2,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of a vector x of dimension dim+1,
 *   with one block of dim threads, on real data.
 *
 * ON ENTRY :
 *   x0hihihi is the highest double of the first element of x; 
 *   x0lohihi is the second highest double of the first element of x; 
 *   x0hilohi is the third highest double of the first element of x; 
 *   x0lolohi is the fourth highest double of the first element of x; 
 *   x0hihilo is the fourth lowest double of the first element of x; 
 *   x0lohilo is the third lowest double of the first element of x; 
 *   x0hilolo is the second lowest double of the first element of x; 
 *   x0lololo is the lowest double of the first element of ; 
 *   x1hihihi is an array with the dim highest doubles of x;
 *   x1lohihi is an array with the dim second highest doubles of x;
 *   x1hilohi is an array with the dim third highest doubles of x;
 *   x1lolohi is an array with the dim fourth highest doubles of x;
 *   x1hihilo is an array with the dim fourth lowest doubles of x;
 *   x1lohilo is an array with the dim third lowest doubles of x;
 *   x1hilolo is an array with the dim second lowest doubles of x;
 *   x1lololo is an array with the dim lowest doubles of x;
 *   dim      the dimension of the vector must equal the block size;
 *   dimLog2  equals ceil(log2((double) dim), used in sum reduction;
 *   vhihihi  has space allocated for dim+1 doubles;
 *   vlohihi  has space allocated for dim+1 doubles;
 *   vhilohi  has space allocated for dim+1 doubles;
 *   vlolohi  has space allocated for dim+1 doubles;
 *   vhihilo  has space allocated for dim+1 doubles;
 *   vlohilo  has space allocated for dim+1 doubles;
 *   vhilolo  has space allocated for dim+1 doubles;
 *   vlololo  has space allocated for dim+1 doubles.
 *
 * ON RETURN :
 *   vhihihi  has the highest doubles of the Householder vector v;
 *   vlohihi  has the second highest doubles of v;
 *   vhilohi  has the third highest doubles of v;
 *   vlolohi  has the fourth highest doubles of v;
 *   vhihilo  has the fourth lowest doubles of v;
 *   vlohilo  has the third lowest doubles of v;
 *   vhilolo  has the second lowest doubles of v;
 *   vlololo  has the lowest doubles of v;
 *   betahihihi is the highest double of 2/(transpose(v)*v);
 *   betalohihi is the second highest double of 2/(transpose(v)*v);
 *   betahilohi is the third highest double of 2/(transpose(v)*v);
 *   betalolohi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v). */

__global__ void cmplx8_small_house
 ( double *x0rehihihi, double *x0relohihi,
   double *x0rehilohi, double *x0relolohi,
   double *x0rehihilo, double *x0relohilo,
   double *x0rehilolo, double *x0relololo,
   double *x0imhihihi, double *x0imlohihi,
   double *x0imhilohi, double *x0imlolohi,
   double *x0imhihilo, double *x0imlohilo,
   double *x0imhilolo, double *x0imlololo,
   double *x1rehihihi, double *x1relohihi,
   double *x1rehilohi, double *x1relolohi,
   double *x1rehihilo, double *x1relohilo,
   double *x1rehilolo, double *x1relololo,
   double *x1imhihihi, double *x1imlohihi,
   double *x1imhilohi, double *x1imlolohi,
   double *x1imhihilo, double *x1imlohilo,
   double *x1imhilolo, double *x1imlololo,
   int dim, int dimLog2,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of a vector x of dimension dim+1,
 *   with one block of dim threads, on complex data.
 *
 * ON ENTRY :
 *   x0rehihihi is the highest double of the real part
 *             of the first element of x; 
 *   x0relohihi is the second highest double of the real part
 *             of the first element of x; 
 *   x0rehilohi is the third highest double of the real part
 *             of the first element of x; 
 *   x0relolohi is the fourth highest double of the real part
 *             of the first element of x; 
 *   x0rehihilo is the fourth lowest double of the real part
 *             of the first element of x; 
 *   x0relohilo is the third lowest double of the real part
 *             of the first element of x; 
 *   x0rehilolo is the second lowest double of the real part
 *             of the first element of x; 
 *   x0relololo is the lowest double of the real part
 *             of the first element of x; 
 *   x0imhihihi is the highest double of the imaginary parts
 *             of the first element of x; 
 *   x0imlohihi is the second highest double of the imaginary parts
 *             of the first element of x; 
 *   x0imhilohi is the third highest double of the imaginary parts
 *             of the first element of x; 
 *   x0imlolohi is the fourth highest double of the imaginary parts
 *             of the first element of x; 
 *   x0imhihilo is the fourth lowest double of the imaginary parts
 *             of the first element of x; 
 *   x0imlohilo is the third lowest double of the imaginary parts
 *             of the first element of x; 
 *   x0imhilolo is the second lowest double of the imaginary parts
 *             of the first element of x; 
 *   x0imlololo is the lowest double of the imaginary parts
 *             of the first element of x; 
 *   x1rehihihi is an array of the dim highest doubles
 *             of the real parts of x;
 *   x1relohihi is an array of the dim second highest doubles
 *             of the real parts of x;
 *   x1rehilohi is an array of the dim third highest doubles
 *             of the real parts of x;
 *   x1relolohi is an array of the dim fourth highest doubles
 *             of the real parts of x;
 *   x1rehihilo is an array of the dim fourth lowest doubles
 *             of the real parts of x;
 *   x1relohilo is an array of the dim third lowest doubles
 *             of the real parts of x;
 *   x1rehilolo is an array of the dim second lowest doubles
 *             of the real parts of x;
 *   x1relololo is an array of the dim lowest doubles
 *             of the real parts of x;
 *   x1imhihihi is an array of the dim highest doubles
 *             of the imaginary parts of x;
 *   x1imlohihi is an array of the dim second highest doubles
 *             of the imaginary parts of x;
 *   x1imhilohi is an array of the dim third highest doubles
 *             of the imaginary parts of x;
 *   x1imlolohi is an array of the dim fourth highest doubles
 *             of the imaginary parts of x;
 *   x1imhihilo is an array of the dim fourth lowest doubles
 *             of the imaginary parts of x;
 *   x1imlohilo is an array of the dim third lowest doubles
 *             of the imaginary parts of x;
 *   x1imhilolo is an array of the dim second lowest doubles
 *             of the imaginary parts of x;
 *   x1imlololo is an array of the dim lowest doubles
 *             of the imaginary parts of x;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction;
 *   vrehihihi has space allocated for dim+1 doubles;
 *   vrelohihi has space allocated for dim+1 doubles;
 *   vrehilohi has space allocated for dim+1 doubles;
 *   vrelolohi has space allocated for dim+1 doubles;
 *   vrehihilo has space allocated for dim+1 doubles;
 *   vrelohilo has space allocated for dim+1 doubles;
 *   vrehilolo has space allocated for dim+1 doubles;
 *   vrelololo has space allocated for dim+1 doubles;
 *   vimhihihi has space allocated for dim+1 doubles;
 *   vimlohihi has space allocated for dim+1 doubles;
 *   vimhilohi has space allocated for dim+1 doubles;
 *   vimlolohi has space allocated for dim+1 doubles;
 *   vimhihilo has space allocated for dim+1 doubles;
 *   vimlohilo has space allocated for dim+1 doubles;
 *   vimhilolo has space allocated for dim+1 doubles;
 *   vimlololo has space allocated for dim+1 doubles.
 *
 * ON RETURN :
 *   vrehihihi has the highest doubles of the real parts
               of the Householder vector v;
 *   vrelohihi has the second highest doubles of the real parts of v;
 *   vrehilohi has the third highest doubles of the real parts of v;
 *   vrelolohi has the fourth highest doubles of the real parts of v;
 *   vrehihilo has the fourth lowest doubles of the real parts of v;
 *   vrelohilo has the third lowest doubles of the real parts of v;
 *   vrehilolo has the second lowest doubles of the real parts of v;
 *   vrelololo has the lowest doubles of the real parts of v;
 *   vimhihihi has the highest doubles of the imaginary parts of v;
 *   vimlohihi has the second highest doubles of the imaginary parts of v;
 *   vimhilohi has the third highest doubles of the imaginary parts of v;
 *   vimlolohi has the fourth highest doubles of the imaginary parts of v;
 *   vimhihilo has the fourth lowest doubles of the imaginary parts of v;
 *   vimlohilo has the third lowest doubles of the imaginary parts of v;
 *   vimhilolo has the second lowest doubles of the imaginary parts of v;
 *   vimlololo has the lowest doubles of the imaginary parts of v;
 *   betahihihi is the highest double of 2/(transpose(v)*v);
 *   betalohihi is the second highest double of 2/(transpose(v)*v);
 *   betahilohi is the third highest double of 2/(transpose(v)*v);
 *   betalohihi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v). */

__global__ void dbl8_large_sum_of_squares
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *sumshihihi, double *sumslohihi,
   double *sumshilohi, double *sumslolohi,
   double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int dim, int BS, int BSLog2 );
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
 *   vhihihi    the highest doubles of a vector v;
 *   vlohihi    the second highest doubles of a vector v;
 *   vhilohi    the third highest doubles of a vector v;
 *   vlolohi    the fourth highest doubles of a vector v;
 *   vhihilo    the fourth lowest doubles of v;
 *   vlohilo    the third lowest doubles of v;
 *   vhilolo    the second lowest doubles of v;
 *   vlololo    the lowest doubles of v;
 *   sumshihihi has space for as many doubles as the number of blocks;
 *   sumslohihi has space for as many doubles as the number of blocks;
 *   sumshilohi has space for as many doubles as the number of blocks;
 *   sumslolohi has space for as many doubles as the number of blocks;
 *   sumshihilo has space for as many doubles as the number of blocks;
 *   sumslohilo has space for as many doubles as the number of blocks;
 *   sumshilolo has space for as many doubles as the number of blocks;
 *   sumslololo has space for as many doubles as the number of blocks;
 *   dim        number of elements in v;
 *   BS         the block size or the number of threads in the block;
 *   BSLog2     equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshihihi are the highest doubles of sums of squares of a vector,
 *              the i-th entry is computed by the i-th block of threads;
 *   sumslohihi are the second highest doubles of the sums of squares;
 *   sumshilohi are the third highest doubles of the sums of squares;
 *   sumslolohi are the fourth highest doubles of the sums of squares;
 *   sumshihilo are the fourth lowest doubles of the sums of squares;
 *   sumslohilo are the third lowest doubles of the sums of squares;
 *   sumshihilo are the second lowest doubles of the sums of squares;
 *   sumslololo are the lowest doubles of the sums of squares. */

__global__ void cmplx8_large_sum_of_squares
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *sumshihihi, double *sumslohihi,
   double *sumshilohi, double *sumslolohi,
   double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int dim, int BS, int BSLog2 );
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
 *   vrehihihi  are the highest doubles of the real parts of v;
 *   vrelohihi  are the second highest doubles of the real parts of v;
 *   vrehilohi  are the third highest doubles of the real parts of v;
 *   vrelolohi  are the fourth highest doubles of the real parts of v;
 *   vrehihilo  are the fourth lowest doubles of the real parts of v;
 *   vrelohilo  are the third lowest doubles of the real parts of v;
 *   vrehilolo  are the second lowest doubles of the real parts of v;
 *   vrelololo  are the lowest doubles of the real parts of v;
 *   vimhihihi  are the highest doubles of the imaginary parts of v;
 *   vimlohihi  are the second highest doubles of the imaginary parts of v;
 *   vimhilohi  are the third highest doubles of the imaginary parts of v;
 *   vimlolohi  are the fourth highest doubles of the imaginary parts of v;
 *   vimhihilo  are the fourth lowest doubles of the imaginary parts of v;
 *   vimlohilo  are the third lowest doubles of the imaginary parts of v;
 *   vimhilolo  are the second lowest doubles of the imaginary parts of v;
 *   vimlololo  are the lowest doubles of the imaginary parts of v;
 *   sumshihihi has space for as many doubles as the number of blocks;
 *   sumslohihi has space for as many doubles as the number of blocks;
 *   sumshilohi has space for as many doubles as the number of blocks;
 *   sumslolohi has space for as many doubles as the number of blocks;
 *   sumshihilo has space for as many doubles as the number of blocks;
 *   sumslohilo has space for as many doubles as the number of blocks;
 *   sumshilolo has space for as many doubles as the number of blocks;
 *   sumslololo has space for as many doubles as the number of blocks;
 *   dim        number of elements in v;
 *   BS         the block size or the number of threads in the block;
 *   BSLog2     equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshihihi are the highest doubles of sums of squares of a vector,
 *              the i-th entry is computed by the i-th block of threads;
 *   sumslohihi are the second highest doubles of the sums of squares;
 *   sumshilohi are the third highest doubles of the sums of squares;
 *   sumslolohi are the fourth highest doubles of the sums of squares;
 *   sumshihilo are the fourth lowest doubles of the sums of squares;
 *   sumslohilo are the third lowest doubles of the sums of squares;
 *   sumshihilo are the second lowest doubles of the sums of squares;
 *   sumslololo are the lowest doubles of the sums of squares. */

__global__ void dbl8_sum_accumulator
 ( double *sumshihihi, double *sumslohihi,
   double *sumshilohi, double *sumslolohi,
   double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo,
   int nbsums, int nbsumsLog2,
   double *acchihihi, double *acclohihi,
   double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo,
   double *acchilolo, double *acclololo );
/*
 * DESCRIPTION :
 *   Accumulates the sum by one block of threads, on real data.
 *
 * REQUIRED : the number of threads should be equal to nbsums.
 *
 * ON ENTRY :
 *   sumshihihi are the highest doubles of sums of squares of vector slices;
 *   sumslohihi are the second highest doubles of sums of squares;
 *   sumshilohi are the third highest doubles of sums of squares;
 *   sumslolohi are the fourth highest doubles of sums of squares;
 *   sumshihilo are the fourth lowest doubles of sums of squares;
 *   sumslohilo are the third lowest doubles of sums of squares;
 *   sumshilolo are the second lowest doubles of sums of squares;
 *   sumslololo are the lowest doubles of sums of squares;
 *   nbsums     the number of elements in sums equals
 *              the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction.
 *
 * ON RETURN :
 *   acchihihi  the highest double of the sum;
 *   acclohihi  the second highest double of the sum;
 *   acchilohi  the third highest double of the sum;
 *   acclolohi  the fourth highest double of the sum;
 *   acchihilo  the fourth lowest double of the sum;
 *   acclohilo  the third lowest double of the sum;
 *   acchilolo  the second lowest double of the sum;
 *   acclololo  the lowest double of the sum. */

__global__ void dbl8_normalize
 ( int dim, int szt,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *v0hihihi, double *v0lohihi, double *v0hilohi, double *v0lolohi,
   double *v0hihilo, double *v0lohilo, double *v0hilolo, double *v0lololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo );
/*
 * DESCRIPTION :
 *   Divides every element in the vector x with the same number v0,
 *   using multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   dim       number of elements in the vectors x and v;
 *   szt       size of one block;
 *   xhihihi   highest doubles of the vector x;
 *   xlohihi   second highest doubles of the vector x;
 *   xhilohi   third highest doubles of the vector x;
 *   xlolohi   fourth highest doubles of the vector x;
 *   xhihilo   fourth lowest doubles of the vector x;
 *   xlohilo   third lowest doubles of the vector x;
 *   xhilolo   second lowest doubles of the vector x;
 *   xlololo   lowest doubles of the vector x;
 *   v0hihihi  highest double of v0;
 *   v0lohihi  second highest double of v0;
 *   v0hilohi  third highest double of v0;
 *   v0lolohi  fourth highest double of v0;
 *   v0hihilo  fourth lowest double of v0;
 *   v0lohilo  third lowest double of v0;
 *   v0hilolo  second lowest double of v0;
 *   v0lololo  lowest double of v0;
 *   vhihihi   space for dim doubles;
 *   vlohihi   space for dim doubles;
 *   vhilohi   space of dim doubles;
 *   vlolohi   space of dim doubles;
 *   vhihilo   space for dim doubles;
 *   vlohilo   space for dim doubles;
 *   vhilolo   space of dim doubles;
 *   vlololo   space of dim doubles.
 *
 * ON RETURN :
 *   vhihihi   highest doubles of x divided by v0;
 *   vlohihi   second highest doubles of x divided by v0;
 *   vhilohi   third highest doubles of x divided by v0;
 *   vlolohi   fourth highest doubles of x divided by v0;
 *   vhihilo   fourth lowest doubles of x divided by v0;
 *   vlohilo   third lowest doubles of x divided by v0;
 *   vhilolo   second lowest doubles of x divided by v0;
 *   vlololo   lowest doubles of x divided by v0. */

__global__ void cmplx8_normalize
 ( int dim, int szt,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *inv0rehihihi, double *inv0relohihi,
   double *inv0rehilohi, double *inv0relolohi,
   double *inv0rehihilo, double *inv0relohilo,
   double *inv0rehilolo, double *inv0relololo,
   double *inv0imhihihi, double *inv0imlohihi,
   double *inv0imhilohi, double *inv0imlolohi,
   double *inv0imhihilo, double *inv0imlohilo,
   double *inv0imhilolo, double *inv0imlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo );
/*
 * DESCRIPTION :
 *   Multiplies every element in the vector x with the same number v0,
 *   using multiple blocks of threads, on complex data.
 *
 * ON ENTRY :
 *   dim       number of elements in the vectors x and v;
 *   szt       size of one block;
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
 *   ximhilohi are the third highest doubles of the imaginary parts of x;
 *   ximlolohi are the fourth highest doubles of the imaginary parts of x;
 *   ximhihilo are the fourth lowest doubles of the imaginary parts of x;
 *   ximlohilo are the third lowest doubles of the imaginary parts of x;
 *   ximhilolo are the second lowest doubles of the imaginary parts of x;
 *   ximlololo are the lowest doubles of the imaginary parts of x;
 *   inv0rehihihi is the highest double of the real part of 1/v0;
 *   inv0relohihi is the second highest double of the real part of 1/v0;
 *   inv0rehilohi is the third highest double of the real part of 1/v0;
 *   inv0relolohi is the fourth highest double of the real part of 1/v0;
 *   inv0rehihilo is the fourth lowest double of the real part of 1/v0;
 *   inv0relohilo is the third lowest double of the real part of 1/v0;
 *   inv0rehilolo is the second lowest double of the real part of 1/v0;
 *   inv0relololo is the lowest double of the real part of 1/v0;
 *   inv0imhihihi is the highest double of the imaginary part of 1/v0;
 *   inv0imlohihi is the second highest double of the imaginary part of 1/v0;
 *   inv0imhilohi is the third highest double of the imaginary part of 1/v0;
 *   inv0imlolohi is the fourth highest double of the imaginary part of 1/v0;
 *   inv0imhihilo is the fourth lowest double of the imaginary part of 1/v0;
 *   inv0imlohilo is the third lowest double of the imaginary part of 1/v0;
 *   inv0imhilolo is the second lowest double of the imaginary part of 1/v0;
 *   inv0imlololo is the lowest double of the imaginary part of 1/v0;
 *   vrehihihi has space for dim doubles;
 *   vrelohihi has space for dim doubles;
 *   vrehilohi has space for dim doubles;
 *   vrelolohi has space for dim doubles;
 *   vrehihilo has space for dim doubles;
 *   vrelohilo has space for dim doubles;
 *   vrehilolo has space for dim doubles;
 *   vrelololo has space for dim doubles;
 *   vimhihihi has space for dim doubles;
 *   vimlohihi has space for dim doubles;
 *   vimhilohi has space for dim doubles;
 *   vimlolohi has space for dim doubles;
 *   vimhihilo has space for dim doubles;
 *   vimlohilo has space for dim doubles;
 *   vimhilolo has space for dim doubles;
 *   vimlololo has space for dim doubles.
 *
 * ON RETURN :
 *   vrehihihi is the highest doubles of x multiplied by 1/v0;
 *   vrelohihi is the second highest doubles of x multiplied by 1/v0;
 *   vrehilohi is the third highest doubles of x multiplied by 1/v0;
 *   vrelolohi is the fourth highest doubles of x multiplied by 1/v0;
 *   vrehihilo is the fourth lowest doubles of x multiplied by 1/v0;
 *   vrelohilo is the third lowest doubles of x multiplied by 1/v0;
 *   vrehilolo is the second lowest doubles of x multiplied by 1/v0;
 *   vrelololo is the lowest doubles of x multiplied by 1/v0;
 *   vimhihihi is the highest doubles of x multiplied by 1/v0;
 *   vimlohihi is the second highest doubles of x multiplied by 1/v0;
 *   vimhilohi is the third highest doubles of x multiplied by 1/v0;
 *   vimlolohi is the fourth highest doubles of x multiplied by 1/v0;
 *   vimhihilo is the fourth lowest doubles of x multiplied by 1/v0;
 *   vimlohilo is the third lowest doubles of x multiplied by 1/v0;
 *   vimhilolo is the second lowest doubles of x multiplied by 1/v0;
 *   vimlololo is the lowest doubles of x multiplied by 1/v0. */

__global__ void cmplx8_normalize_rere
 ( int dim, int szt,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *inv0rehihihi, double *inv0relohihi,
   double *inv0rehilohi, double *inv0relolohi,
   double *inv0rehihilo, double *inv0relohilo,
   double *inv0rehilolo, double *inv0relololo,
   double *vrehihihi, double *vrelohihi,
   double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo,
   double *vrehilolo, double *vrelololo );
/*
 * DESCRIPTION :
 *   Does the first multiplication in the normalization. */

__global__ void cmplx8_normalize_imim
 ( int dim, int szt,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *inv0imhihihi, double *inv0imlohihi,
   double *inv0imhilohi, double *inv0imlolohi,
   double *inv0imhihilo, double *inv0imlohilo,
   double *inv0imhilolo, double *inv0imlololo,
   double *vrehihihi, double *vrelohihi,
   double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo,
   double *vrehilolo, double *vrelololo );
/*
 * DESCRIPTION :
 *   Does the second multiplication in the normalization. */

__global__ void cmplx8_normalize_imre
 ( int dim, int szt,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *inv0rehihihi, double *inv0relohihi,
   double *inv0rehilohi, double *inv0relolohi,
   double *inv0rehihilo, double *inv0relohilo,
   double *inv0rehilolo, double *inv0relololo,
   double *vimhihihi, double *vimlohihi,
   double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo );
/*
 * DESCRIPTION :
 *   Does the third multiplication in the normalization. */

__global__ void cmplx8_normalize_reim
 ( int dim, int szt,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *inv0imhihihi, double *inv0imlohihi,
   double *inv0imhilohi, double *inv0imlolohi,
   double *inv0imhihilo, double *inv0imlohilo,
   double *inv0imhilolo, double *inv0imlololo,
   double *vimhihihi, double *vimlohihi,
   double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo );
/*
 * DESCRIPTION :
 *   Does the fourth multiplication in the normalization. */

__global__ void dbl8_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo );
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
 *   Rhihihi  highest doubles of R, stored column wise;
 *   Rlohihi  second highest doubles of R, stored column wise;
 *   Rhilohi  third highest doubles of R, stored column wise;
 *   Rlolohi  fourth highest doubles of R, stored column wise;
 *   Rhihilo  fourth lowest doubles of R, stored column wise;
 *   Rlohilo  third lowest doubles of R, stored column wise;
 *   Rhihilo  second lowest doubles of R, stored column wise;
 *   Rlololo  lowest doubles of R, stored column wise;
 *   vhihihi  highest doubles of the Householder vector;
 *   vlohihi  second highest doubles of the Householder vector;
 *   vhilohi  third highest doubles of the Householder vector;
 *   vlolohi  fourth highest doubles of the Householder vector;
 *   vhihilo  fourth lowest doubles of the Householder vector;
 *   vlohilo  third lowest doubles of the Householder vector;
 *   vhilolo  second lowest doubles of the Householder vector;
 *   vlololo  lowest doubles of the Householder vector;
 *   betalohihi is the second highest double of 2/(transpose(v)*v);
 *   betahilohi is the third highest double of 2/(transpose(v)*v);
 *   betalohihi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v).
 *
 * ON RETURN :
 *   Rhihihi  highest doubles of the updated R, which is trapezoidal;
 *   Rlohihi  second highest doubles of the updated R;
 *   Rhilohi  third highest doubles of the updated R;
 *   Rlolohi  fourth highest doubles of the updated R;
 *   Rhihilo  fourth lowest doubles of the updated R;
 *   Rlohilo  third lowest doubles of the updated R;
 *   Rhilolo  second lowest doubles of the updated R;
 *   Rlololo  lowest doubles of the updated R. */

__global__ void cmplx8_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi, 
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo );
/*
 * DESCRIPTION :
 *   Updates the matrix R starting at column k
 *   with the Householder vector in v and beta, on complex data,
 *   with a total of ncols - k threads in one block.
 *
 * ON ENTRY :
 *   nrows     number of rows of R;
 *   ncols     number of columns of R;
 *   szt       size of each block;
 *   k         index of the current column;
 *   Rrehihihi are the highest doubles of the real parts of R,
 *             stored column wise;
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
 *   Rimlohihi are the fourth highest doubles of the imaginary parts of R;
 *   Rimhihilo are the fourth lowest doubles of the imaginary parts of R;
 *   Rimlohilo are the third lowest doubles of the imaginary parts of R;
 *   Rimhilolo are the second lowest doubles of the imaginary parts of R;
 *   Rimlololo are the lowest doubles of the imaginary parts of R;
 *   betalohihi is the second highest double of 2/(transpose(v)*v);
 *   betahilohi is the third highest double of 2/(transpose(v)*v);
 *   betalohihi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v).
 *
 * ON RETURN :
 *   Rrehihihi are the highest doubles of the real parts of the updated R;
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
 *   Rimlololo are the lowest doubles of the imaginary parts of R. */

__global__ void dbl8_RTdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *RTdotvhihihi, double *RTdotvlohihi,
   double *RTdotvhilohi, double *RTdotvlolohi,
   double *RTdotvhihilo, double *RTdotvlohilo,
   double *RTdotvhilolo, double *RTdotvlololo );
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
 *   Rhihihi  highest doubles of a column wise stored matrix R
 *            with number of rows equal to nrows;
 *   Rlohihi  second highest doubles of R;
 *   Rhilohi  third highest doubles of R;
 *   Rlolohi  fourth highest doubles of R;
 *   Rhihilo  fourth lowest doubles of R;
 *   Rlohilo  third lowest doubles of R;
 *   Rhilolo  second lowest doubles of R;
 *   Rlololo  lowest doubles of R;
 *   vhihi    start of the highest doubles of the first nonzero element
 *            of a Householder vector v;
 *   vlohihi  start of the second highest doubles of v;
 *   vhilohi  start of the third highest doubles of v;
 *   vlolohi  start of the fourth highest doubles of v;
 *   vhihilo  start of the fourth lowest doubles of v;
 *   vlohilo  start of the third lowest doubles of v;
 *   vhilolo  start of the second lowest doubles of v;
 *   vlololo  start of the lowest doubles of v;
 *   RTdotvhihihi has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvlohihi has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvhilohi has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvlolohi has space for a matrix of nrows-by-szt, plus some padding.
 *   RTdotvhihilo has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvlohilo has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvhilolo has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvlololo has space for a matrix of nrows-by-szt, plus some padding.
 *
 * ON RETURN :
 *   RTdotvhihihi are the highest doubles of the element-by-element products
 *             of R^T with v, stored row by row;
 *   RTdotvlohihi are the second highest doubles of RTdotv;
 *   RTdotvhilohi are the third highest doubles of RTdotv;
 *   RTdotvlolohi are the fourth highest doubles of RTdotv;
 *   RTdotvhihilo are the fourth lowest doubles of RTdotv;
 *   RTdotvlohilo are the third lowest doubles of RTdotv;
 *   RTdotvhihilo are the second lowest doubles of RTdotv;
 *   RTdotvlololo are the lowest doubles of RTdotv. */

__global__ void cmplx8_RHdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo );
/*
 * DESCRIPTION :
 *   The elements of the matrix RHdotv are the elements of R^H,
 *   multiplied with the corresponding element of v.
 *   Multiple blocks of threads operate on complex data.
 *
 * ON ENTRY :
 *   nrows     number of rows of R;
 *   szt       size of one tile and number of threads in one block;
 *   colidx    index of the current column in R;
 *   Roffset   offset in R for the first row to start;
 *   dim       number of columns in R;
 *   Rrehihihi are the highest doubles of the real parts of R, a column wise
               stored matrix with number of rows equal to nrows;
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

__global__ void cmplx8_RHdotvRe
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo );
/*
 * DESCRIPTION :
 *   Computes only the real parts of RHdotv. */

__global__ void cmplx8_RHdotvReRe
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo );
/*
 * DESCRIPTION :
 *   Does the multiplication of the real parts only to update RHdotvre. */

__global__ void cmplx8_RHdotvImIm
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo );
/*
 * DESCRIPTION :
 *   Updates RHdotvre with the multiplication of the imaginary parts. */

__global__ void cmplx8_RHdotvIm
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo );
/*
 * DESCRIPTION :
 *   Computes only the imaginary parts of RHdotv. */

__global__ void cmplx8_RHdotvReIm
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo );
/*
 * DESCRIPTION :
 *   Multiplies only the real parts of R with imaginary parts of v. */

__global__ void cmplx8_RHdotvImRe
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo );
/*
 * DESCRIPTION :
 *   Updates RHdotvim with the imaginary parts of R multiplied
 *   with the real parts of v. */

__global__ void dbl8_sum_betaRTdotv
 ( int nrows,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *RTdotvhihihi, double *RTdotvlohihi,
   double *RTdotvhilohi, double *RTdotvlolohi,
   double *RTdotvhihilo, double *RTdotvlohilo,
   double *RTdotvhilolo, double *RTdotvlololo,
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo );
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

__global__ void cmplx8_sum_betaRHdotv
 ( int nrows,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *RTdotvrehihihi, double *RTdotvrelohihi,
   double *RTdotvrehilohi, double *RTdotvrelolohi,
   double *RTdotvrehihilo, double *RTdotvrelohilo,
   double *RTdotvrehilolo, double *RTdotvrelololo,
   double *RTdotvimhihihi, double *RTdotvimlohihi,
   double *RTdotvimhilohi, double *RTdotvimlolohi,
   double *RTdotvimhihilo, double *RTdotvimlohilo,
   double *RTdotvimhilolo, double *RTdotvimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo );
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

__global__ void dbl8_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo );
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

__global__ void cmplx8_medium_subvbetaRHv
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo );
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

__global__ void cmplx8_medium_subvbetaRHvRe
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo );
/*
 * DESCRIPTION :
 *   Updates only the real parts of R. */

__global__ void cmplx8_medium_subvbetaRHvReRe
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo );
/*
 * DESCRIPTION :
 *   Multiplies only the real parts of R and v to update R. */

__global__ void cmplx8_medium_subvbetaRHvImIm
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *wimhihihi, double *wimlohihi, double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo, double *wimhilolo, double *wimlololo );
/*
 * DESCRIPTION :
 *   Multiplies only the imaginary parts of R and v to update R. */

__global__ void cmplx8_medium_subvbetaRHvIm
 ( int nrows, int ncols, int szt, int k,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo );
/*
 * DESCRIPTION :
 *   Updates only the imaginary parts of R. */

__global__ void cmplx8_medium_subvbetaRHvImRe
 ( int nrows, int ncols, int szt, int k,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo );
/*
 * DESCRIPTION :
 *   Updates the imaginary parts of R by the multiplication of imaginary
 *   parts with real parts. */

__global__ void cmplx8_medium_subvbetaRHvReIm
 ( int nrows, int ncols, int szt, int k,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *wimhihihi, double *wimlohihi, double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo, double *wimhilolo, double *wimlololo );
/*
 * DESCRIPTION :
 *   Updates the imaginary parts of R by the multiplication of the real
 *   parts with imaginary parts. */

__global__ void dbl8_beta_times_V
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo );
/*
 * DESCRIPTION :
 *   Computes the first vector in the W representation of the Householder
 *   transformations, multiplying B[0] with the first vector of V,
 *   and flipping the sign, with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V and W;
 *   szt      size of one tile and the number of threads in a block;
 *   Bhihihi  highest double of the betas computed by house;
 *   Blohihi  second highest double of the betas;
 *   Bhilohi  third highest double of the betas;
 *   Blolohi  fourth highest double of the betas;
 *   Bhihilo  fourth lowest double of the betas;
 *   Blohilo  third lowest double of the betas;
 *   Bhilolo  second lowest double of the betas;
 *   Blololo  lowest double of the betas;
 *   Vhihihi  Vhihihi[nrows*i] is the first highest doubles of the i-th
 *            Householder vector, with i zeros inserted so V is trapezoidal;
 *   Vlohihi  second highest doubles of V;
 *   Vhilohi  third highest doubles of V;
 *   Vlolohi  fourth highest doubles of V;
 *   Vhihilo  fourth lowest doubles of V;
 *   Vlohilo  third lowest doubles of V;
 *   Vhilolo  second lowest doubles of V;
 *   Vlololo  lowest doubles of V;
 *   Whihihi  space for the highest doubles of W;
 *   Wlohihi  space for the second highest doubles of W;
 *   Whilohi  space for the third highest doubles of W;
 *   Wlolohi  space for the fourth highest doubles of W;
 *   Whihilo  space for the fourth lowest doubles of W;
 *   Wlohilo  space for the third lowest doubles of W;
 *   Whilolo  space for the second lowest doubles of W;
 *   Wlololo  space for the lowest doubles of W.
 *
 * ON RETURN :
 *   Whihihi  the first nrows numbers store the highest doubles of the
 *            first vector of the W matrix in the WY representation;
 *   Wlohihi  the second highest doubles of W;
 *   Whilohi  the third highest doubles of W;
 *   Wlolohi  the fourth highest doubles of W;
 *   Whihilo  the fourth lowest doubles of W;
 *   Wlohilo  the third lowest doubles of W;
 *   Whilolo  the second lowest doubles of W;
 *   Wlololo  the lowest doubles of W. */

__global__ void cmplx8_beta_times_V
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi,
   double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo,
   double *Wimhilolo, double *Wimlololo );
/*
 * DESCRIPTION :
 *   Computes the first vector in the W representation of the Householder
 *   transformations, multiplying B[0] with the first vector of V,
 *   and flipping the sign, with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   nrows     number of rows in V and W;
 *   szt       size of one tile and the number of threads in a block;
 *   Bhihihi   highest double of the betas computed by house;
 *   Blohihi   second highest double of the betas;
 *   Bhilohi   third highest double of the betas;
 *   Blolohi   fourth highest double of the betas;
 *   Bhihilo   fourth lowest double of the betas;
 *   Blohilo   third lowest double of the betas;
 *   Bhilolo   second lowest double of the betas;
 *   Blololo   lowest double of the betas;
 *   Vrehihihi Vrehihihi[nrows*i] is the first highest doubles of the real
 *             parts of the i-th Householder vector, with i zeros inserted;
 *   Vrelohihi are the second highest doubles of the real parts of V;
 *   Vrehilohi are the third highest doubles of the real parts of V;
 *   Vrelolohi are the fourth highest doubles of the real parts of V;
 *   Vrehihilo are the fourth lowest doubles of the real parts of V;
 *   Vrelohilo are the third lowest doubles of the real parts of V;
 *   Vrehilolo are the second lowest doubles of the real parts of V;
 *   Vrelololo are the lowest doubles of the real parts of V;
 *   Vimhihihi are the highest doubles of the imaginary parts of V;
 *   Vimlohihi are the second highest doubles of the imaginary parts of V;
 *   Vimhilohi are the third highest doubles of the imaginary parts of V;
 *   Vimlolohi are the fourth highest doubles of the imaginary parts of V;
 *   Vimhihilo are the fourth lowest doubles of the imaginary parts of V;
 *   Vimlohilo are the third lowest doubles of the imaginary parts of V;
 *   Vimhilolo are the second lowest doubles of the imaginary parts of V;
 *   Vimlololo are the lowest doubles of the imaginary parts of V;
 *   Wrehihihi has space for the highest doubles of the real parts of the W;
 *   Wrelohihi has space for the second highest doubles
 *             of the real parts of the W;
 *   Wrehilohi has space for the third highest doubles
 *             of the real parts of the W;
 *   Wrelolohi has space for the fourth highest doubles
 *             of the real parts of the W;
 *   Wrehihilo has space for the fourth lowest doubles
 *             of the real parts of the W;
 *   Wrelohilo has space for the third lowest doubles
 *             of the real parts of the W;
 *   Wrehilolo has space for the second lowest doubles
 *             of the real parts of the W;
 *   Wrelololo has space for the lowest doubles of the real parts of the W;
 *   Wimhihihi has space for the highest doubles of the imaginary parts of W;
 *   Wimlohihi has space for the second highest doubles
 *             of the imaginary parts of W;
 *   Wimhilohi has space for the third highest doubles
 *             of the imaginary parts of W;
 *   Wimlolohi has space for the fourth highest doubles
 *             of the imaginary parts of W;
 *   Wimhihilo has space for the fourth lowest doubles
 *             of the imaginary parts of W.
 *   Wimlohilo has space for the third lowest doubles
 *             of the imaginary parts of W.
 *   Wimhilolo has space for the second lowest doubles
 *             of the imaginary parts of W.
 *   Wimlololo has space for the lowest doubles of the imaginary parts of W.
 *
 * ON RETURN :
 *   Wrehihihi has the highest doubles of the real in the first nrows numbers
 *             parts of the first vector of the W the WY representation;
 *   Wrelohihi has the second highest doubles of the real parts of W;
 *   Wrehilohi has the third highest doubles of the real parts of W;
 *   Wrelolohi has the fourth highest doubles of the real parts of W;
 *   Wrehihilo has the fourth lowest doubles of the real parts of W;
 *   Wrelohilo has the third lowest doubles of the real parts of W;
 *   Wrehilolo has the second lowest doubles of the real parts of W;
 *   Wrelololo has the lowest doubles of the real parts of W;
 *   Wimhihihi has the highest doubles of the imaginary parts of W;
 *   Wimlohihi has the second highest doubles of the imaginary parts of W;
 *   Wimhilohi has the third highest doubles of the imaginary parts of W;
 *   Wimlolohi has the fourth highest doubles of the imaginary parts of W;
 *   Wimhihilo has the fourth lowest doubles of the imaginary parts of W;
 *   Wimlohilo has the third lowest doubles of the imaginary parts of W;
 *   Wimhilolo has the second lowest doubles of the imaginary parts of W;
 *   Wimlololo has the lowest doubles of the imaginary parts of W. */

__global__ void dbl8_initialize_WYT
 ( int dim, int szt,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo );
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

__global__ void cmplx8_initialize_WYH
 ( int dim, int szt,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo );
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

__global__ void dbl8_update_WYT
 ( int dim, int szt,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo );
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

__global__ void cmplx8_update_WYH
 ( int dim, int szt,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYHrehihihi, double *WYHrelohihi,
   double *WYHrehilohi, double *WYHrelolohi,
   double *WYHrehihilo, double *WYHrelohilo,
   double *WYHrehilolo, double *WYHrelololo,
   double *WYHimhihihi, double *WYHimlohihi,
   double *WYHimhilohi, double *WYHimlolohi,
   double *WYHimhihilo, double *WYHimlohilo,
   double *WYHimhilolo, double *WYHimlololo );
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

__global__ void dbl8_beta_next_W
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo );
/*
 * DECRIPTION :
 *   Computes the next column in the W matrix, with multiple blocks,
 *   on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V, W, and the dimension of WYT;
 *   szt      number of threads in one block;
 *   Bhi      highest double of the beta with the next Householder vector;
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

__global__ void cmplx8_beta_next_W
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYHrehihihi, double *WYHrelohihi,
   double *WYHrehilohi, double *WYHrelolohi,
   double *WYHrehihilo, double *WYHrelohilo,
   double *WYHrehilolo, double *WYHrelololo,
   double *WYHimhihihi, double *WYHimlohihi,
   double *WYHimhilohi, double *WYHimlolohi,
   double *WYHimhihilo, double *WYHimlohilo,
   double *WYHimhilolo, double *WYHimlololo );
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

__global__ void dbl8_small_WYT
 ( int nrows, int szt,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo );
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

__global__ void cmplx8_small_WYH
 ( int nrows, int szt,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *Yrehihihi, double *Yrelohihi, double *Yrehilohi, double *Yrelolohi,
   double *Yrehihilo, double *Yrelohilo, double *Yrehilolo, double *Yrelololo,
   double *Yimhihihi, double *Yimlohihi, double *Yimhilohi, double *Yimlolohi,
   double *Yimhihilo, double *Yimlohilo, double *Yimhilolo, double *Yimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo );
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

__global__ void dbl8_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihihi, double *Qlohihi, double *Qhilohi, double *Qlolohi,
   double *Qhihilo, double *Qlohilo, double *Qhilolo, double *Qlololo,
   double *WYThihihi, double *WYTlohihi, double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo, double *WYThilolo, double *WYTlololo,
   double *QWYThihihi, double *QWYTlohihi,
   double *QWYThilohi, double *QWYTlolohi,
   double *QWYThihilo, double *QWYTlohilo,
   double *QWYThilolo, double *QWYTlololo );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYT into the matrix QWYT, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYT matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qhihihi  the highest doubles of Q;
 *   Qlohihi  the second highest doubles of Q;
 *   Qhilohi  the third highest doubles of Q;
 *   Qlolohi  the fourth highest doubles of Q;
 *   Qhihilo  the fourth lowest doubles of Q;
 *   Qlohilo  the third lowest doubles of Q;
 *   Qhilolo  the second lowest doubles of Q;
 *   Qlololo  the lowest doubles of Q;
 *   WYThihihi are the highest doubles of W*Y^T;
 *   WYTlohihi are the second highest doubles of W*Y^T;
 *   WYThilohi are the third highest doubles of W*Y^T;
 *   WYTlohihi are the fourth highest doubles of W*Y^T;
 *   WYThihilo are the fourth lowest doubles of W*Y^T;
 *   WYTlohilo are the third lowest doubles of W*Y^T;
 *   WYThilolo are the second lowest doubles of W*Y^T;
 *   WYTlololo are the lowest doubles of W*Y^T;
 *   QWYThihihi has space for a dim-by-dim matrix;
 *   QWYTlohihi has space for a dim-by-dim matrix;
 *   QWYThilohi has space for a dim-by-dim matrix;
 *   QWYTlolohi has space for a dim-by-dim matrix;
 *   QWYThihilo has space for a dim-by-dim matrix;
 *   QWYTlohilo has space for a dim-by-dim matrix;
 *   QWYThilolo has space for a dim-by-dim matrix;
 *   QWYTlololo has space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYThihihi are the highest doubles of Q*WYT;
 *   QWYTlohihi are the second highest doubles of Q*WYT;
 *   QWYThilohi are the third highest doubles of Q*WYT;
 *   QWYTlolohi are the fourth highest doubles of Q*WYT;
 *   QWYThihilo are the fourth lowest doubles of Q*WYT;
 *   QWYTlohilo are the third lowest doubles of Q*WYT;
 *   QWYThilolo are the second lowest doubles of Q*WYT;
 *   QWYTlololo are the lowest doubles of Q*WYT. */

__global__ void cmplx8_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihihi, double *Qrelohihi, double *Qrehilohi, double *Qrelolohi,
   double *Qrehihilo, double *Qrelohilo, double *Qrehilolo, double *Qrelololo,
   double *Qimhihihi, double *Qimlohihi, double *Qimhilohi, double *Qimlolohi,
   double *Qimhihilo, double *Qimlohilo, double *Qimhilolo, double *Qimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo,
   double *QWYTrehihihi, double *QWYTrelohihi,
   double *QWYTrehilohi, double *QWYTrelolohi,
   double *QWYTrehihilo, double *QWYTrelohilo,
   double *QWYTrehilolo, double *QWYTrelololo,
   double *QWYTimhihihi, double *QWYTimlohihi,
   double *QWYTimhilohi, double *QWYTimlolohi,
   double *QWYTimhihilo, double *QWYTimlohilo,
   double *QWYTimhilolo, double *QWYTimlololo );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYT into the matrix QWYT, on complex data.
 *
 * ON ENTRY :
 *   dim       number of rows and columns of the Q matrix;
 *   rowdim    number of rows and columns of the WYT matrix;
 *   szt       the number of threads in a block;
 *   coloff    offset for the column index in QWYT;
 *   Qrehihihi are the highest doubles of the real parts of Q;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third highest doubles of the real parts of Q;
 *   Qrelolohi are the fourth highest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest double of the real parts of Q;
 *   Qrelohilo are the third lowest double of the real parts of Q;
 *   Qrehilolo are the second lowest double of the real parts of Q;
 *   Qrelololo are the lowest double of the real parts of Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi are the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q;
 *   WYTrehihihi are the highest doubles of the real parts of W*Y^T;
 *   WYTrelohihi are the second highest doubles of the real parts of W*Y^T;
 *   WYTrehilohi are the third highest doubles of the real parts of W*Y^T;
 *   WYTrelolohi are the fourth highest doubles of the real parts of W*Y^T;
 *   WYTrehihilo are the fourth lowest doubles of the real parts of W*Y^T;
 *   WYTrelohilo are the third lowest doubles of the real parts of W*Y^T;
 *   WYTrehilolo are the second lowest doubles of the real parts of W*Y^T;
 *   WYTrelololo are the lowest doubles of the real parts of W*Y^T;
 *   WYTimhihihi are the highest doubles of the imaginary parts of W*Y^T;
 *   WYTimlohihi are the second highest doubles of the imaginary parts of W*Y^T;
 *   WYTimhilohi are the third highest doubles of the imaginary parts of W*Y^T;
 *   WYTimlolohi are the fourth highest doubles of the imaginary parts of W*Y^T;
 *   WYTimhihilo are the fourth lowest doubles of the imaginary parts of W*Y^T;
 *   WYTimlohilo are the third lowest doubles of the imaginary parts of W*Y^T;
 *   WYTimhilolo are the second lowest doubles of the imaginary parts of W*Y^T;
 *   WYTimlololo are the lowest doubles of the imaginary parts of W*Y^T;
 *   QWYTrehihihi has space for a dim-by-dim matrix;
 *   QWYTrelohihi has space for a dim-by-dim matrix;
 *   QWYTrehilohi has space for a dim-by-dim matrix;
 *   QWYTrelolohi has space for a dim-by-dim matrix;
 *   QWYTrehihilo has space for a dim-by-dim matrix;
 *   QWYTrelohilo has space for a dim-by-dim matrix;
 *   QWYTrehilolo has space for a dim-by-dim matrix;
 *   QWYTrelololo has space for a dim-by-dim matrix;
 *   QWYTimhihihi has space for a dim-by-dim matrix;
 *   QWYTimlohihi has space for a dim-by-dim matrix;
 *   QWYTimhilohi has space for a dim-by-dim matrix;
 *   QWYTimlolohi has space for a dim-by-dim matrix;
 *   QWYTimhihilo has space for a dim-by-dim matrix;
 *   QWYTimlohilo has space for a dim-by-dim matrix;
 *   QWYTimhilolo has space for a dim-by-dim matrix;
 *   QWYTimlololo has space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYTrehihihi are the highest doubles of the real parts of Q*W*Y^T;
 *   QWYTrelohihi are the second highest doubles of the real parts of Q*W*Y^T;
 *   QWYTrehilohi are the third highest doubles of the real parts of Q*W*Y^T;
 *   QWYTrelolohi are the fourth highest doubles of the real parts of Q*W*Y^T;
 *   QWYTrehihilo are the fourth lowest doubles of the real parts of Q*W*Y^T;
 *   QWYTrelohilo are the third lowest doubles of the real parts of Q*W*Y^T;
 *   QWYTrehilolo are the second lowest doubles of the real parts of Q*W*Y^T;
 *   QWYTrelololo are the lowest doubles of the real parts of Q*W*Y^T;
 *   QWYTimhihihi are the second highest doubles of the imaginary parts
 *                of Q*W*Y^T;
 *   QWYTimlohihi are the third highest doubles of the imaginary parts
 *                of Q*W*Y^T;
 *   QWYTimlolohi are the fourth highest doubles of the imaginary parts 
 *                of Q*W*Y^T;
 *   QWYTimhilolo are the fourth lowest doubles of the imaginary parts
 *                of Q*W*Y^T;
 *   QWYTimlohilo are the third lowest doubles of the imaginary parts
 *                of Q*W*Y^T;
 *   QWYTimhilolo are the second lowest doubles of the imaginary parts
 *                of Q*W*Y^T;
 *   QWYTimlololo are the lowest doubles of the imaginary parts
 *                of Q*W*Y^T. */

__global__ void dbl8_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWThihihi, double *YWTlohihi, double *YWThilohi, double *YWTlolohi,
   double *YWThihilo, double *YWTlohilo, double *YWThilolo, double *YWTlololo,
   double *Chihihi, double *Clohihi, double *Chilohi, double *Clolohi,
   double *Chihilo, double *Clohilo, double *Chilolo, double *Clololo,
   double *YWTChihihi, double *YWTClohihi,
   double *YWTChilohi, double *YWTClolohi,
   double *YWTChihilo, double *YWTClohilo,
   double *YWTChilolo, double *YWTClololo );
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

__global__ void cmplx8_small_YWHC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWTrehihihi, double *YWTrelohihi,
   double *YWTrehilohi, double *YWTrelolohi,
   double *YWTrehihilo, double *YWTrelohilo,
   double *YWTrehilolo, double *YWTrelololo,
   double *YWTimhihihi, double *YWTimlohihi,
   double *YWTimhilohi, double *YWTimlolohi,
   double *YWTimhihilo, double *YWTimlohilo,
   double *YWTimhilolo, double *YWTimlololo,
   double *Crehihihi, double *Crelohihi, double *Crehilohi, double *Crelolohi,
   double *Crehihilo, double *Crelohilo, double *Crehilolo, double *Crelololo,
   double *Cimhihihi, double *Cimlohihi, double *Cimhilohi, double *Cimlolohi,
   double *Cimhihilo, double *Cimlohilo, double *Cimhilolo, double *Cimlololo,
   double *YWTCrehihihi, double *YWTCrelohihi,
   double *YWTCrehilohi, double *YWTCrelolohi,
   double *YWTCrehihilo, double *YWTCrelohilo,
   double *YWTCrehilolo, double *YWTCrelololo,
   double *YWTCimhihihi, double *YWTCimlohihi,
   double *YWTCimhilohi, double *YWTCimlolohi,
   double *YWTCimhihilo, double *YWTCimlohilo,
   double *YWTCimhilolo, double *YWTCimlololo );
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
 *   YWTCimlohi has space for a rowdim-by-coldim matrix.
 *   YWTCimhilo has space for a rowdim-by-coldim matrix.
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

__global__ void dbl8_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihihi, double *Qlohihi, double *Qhilohi, double *Qlolohi,
   double *Qhihilo, double *Qlohilo, double *Qhilolo, double *Qlololo,
   double *QWYThihihi, double *QWYTlohihi,
   double *QWYThilohi, double *QWYTlolohi,
   double *QWYThihilo, double *QWYTlohilo,
   double *QWYThilolo, double *QWYTlololo );
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

__global__ void cmplx8_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihihi, double *Qrelohihi, double *Qrehilohi, double *Qrelolohi,
   double *Qrehihilo, double *Qrelohilo, double *Qrehilolo, double *Qrelololo,
   double *Qimhihihi, double *Qimlohihi, double *Qimhilohi, double *Qimlolohi,
   double *Qimhihilo, double *Qimlohilo, double *Qimhilolo, double *Qimlololo,
   double *QWYTrehihihi, double *QWYTrelohihi,
   double *QWYTrehilohi, double *QWYTrelolohi,
   double *QWYTrehihilo, double *QWYTrelohilo,
   double *QWYTrehilolo, double *QWYTrelololo,
   double *QWYTimhihihi, double *QWYTimlohihi,
   double *QWYTimhilohi, double *QWYTimlolohi,
   double *QWYTimhihilo, double *QWYTimlohilo,
   double *QWYTimhilolo, double *QWYTimlololo );
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

__global__ void dbl8_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *YWTChihihi, double *YWTClohihi,
   double *YWTChilohi, double *YWTClolohi,
   double *YWTChihilo, double *YWTClohilo,
   double *YWTChilolo, double *YWTClololo );
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
 *   YWTChilo are the second lowest doubles of YWT*C.
 *   YWTClolo are the lowest doubles of YWT*C.
 *
 * ON RETURN :
 *   Rhihi    highest doubles of R + YWTC;
 *   Rlohi    second highest doubles of R + YWTC;
 *   Rhilo    second lowest doubles of R + YWTC;
 *   Rlolo    lowest doubles of R + YWTC. */

__global__ void cmplx8_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *YWTCrehihihi, double *YWTCrelohihi,
   double *YWTCrehilohi, double *YWTCrelolohi,
   double *YWTCrehihilo, double *YWTCrelohilo,
   double *YWTCrehilolo, double *YWTCrelololo,
   double *YWTCimhihihi, double *YWTCimlohihi,
   double *YWTCimhilohi, double *YWTCimlolohi,
   double *YWTCimhihilo, double *YWTCimlohilo,
   double *YWTCimhilolo, double *YWTCimlololo );
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

void GPU_dbl8_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *vhihihi_h, double *vlohihi_h, double *vhilohi_h, double *vlolohi_h,
   double *vhihilo_h, double *vlohilo_h, double *vhilolo_h, double *vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the Householder vector for small
 *   enough matrices for the vector to fit entirely in shared memory.
 *
 * ON ENTRY :
 *   nrows     number of rows in the matrix A;
 *   ncols     number of columns in the matrix A;
 *   szt       size of one tile;
 *   nbt       number of tiles, szt*nbt = ncols;
 *   colidx    global index of the current column;
 *   nrows1    number of threads in the block equals the number
 *             of elements computed in the Householder vector;
 *   L         local index of the column in the current tile;
 *   Ahihihi_h are the highest doubles of the matrix on the host;
 *   Alohihi_h are the second highest doubles of the matrix on the host;
 *   Ahilohi_h are the third highest doubles of the matrix on the host;
 *   Alolohi_h are the fourth highest doubles of the matrix on the host;
 *   Ahihilo_h are the fourth lowest doubles of the matrix on the host;
 *   Alohilo_h are the third lowest doubles of the matrix on the host;
 *   Ahilolo_h are the second lowest doubles of the matrix on the host;
 *   Alololo_h are the lowest doubles of the matrix on the host;
 *   Ahihihi_d are the highest doubles of the matrix on the device;
 *   Alohihi_d are the second highest doubles of the matrix on the device;
 *   Ahilohi_d are the third highest doubles of the matrix on the device;
 *   Alolohi_d are the fourth highest doubles of the matrix on the device;
 *   Ahihilo_d are the fourth lowest doubles of the matrix on the device;
 *   Alohilo_d are the third lowest doubles of the matrix on the device;
 *   Ahilolo_d are the second lowest doubles of the matrix on the device;
 *   Alololo_d are the lowest doubles of the matrix on the device;
 *   vhihihi_h has space for the current Householder vector;
 *   vlohihi_h has space for the current Householder vector;
 *   vhilohi_h has space for the current Householder vector;
 *   vlolohi_h has space for the current Householder vector;
 *   vhihilo_h has space for the current Householder vector;
 *   vlohilo_h has space for the current Householder vector;
 *   vhilolo_h has space for the current Householder vector;
 *   vlololo_h has space for the current Householder vector;
 *   Vhihihi_d has space for the Householder vectors on the device;
 *   Vlohihi_d has space for the Householder vectors on the device;
 *   Vhilohi_d has space for the Householder vectors on the device;
 *   Vlolohi_d has space for the Householder vectors on the device;
 *   Vhihilo_d has space for the Householder vectors on the device;
 *   Vlohilo_d has space for the Householder vectors on the device;
 *   Vhilolo_d has space for the Householder vectors on the device;
 *   Vlololo_d has space for the Householder vectors on the device;
 *   betahihihi_h has space for the highest doubles of the betas, if verbose;
 *   betalohihi_h has space for the second highest doubles of the betas,
 *             if verbose;
 *   betahilohi_h has space for the third highest doubles of the betas,
 *             if verbose;
 *   betalolohi_h has space for the fourth highest doubles of the betas,
 *             if verbose;
 *   betahihilo_h has space for the fourth lowest doubles of the betas,
 *             if verbose;
 *   betalohilo_h has space for the third lowest doubles of the betas,
 *             if verbose;
 *   betahilolo_h has space for the second lowest doubles of the betas,
 *             if verbose;
 *   betalololo_h has space for the lowest doubles of the betas, if verbose;
 *   betahihihi_d has space on the device for the highest doubles
 *             of the betas;
 *   betalohihi_d has space on the device for the second highest doubles
 *             of the betas;
 *   betahilohi_d has space on the device for the third highest doubles
 *             of the betas;
 *   betalolohi_d has space on the device for the fourth highest doubles
 *             of the betas;
 *   betahihilo_d has space on the device for the fourth lowest doubles
 *             of the betas;
 *   betalohilo_d has space on the device for the third lowest doubles
 *             of the betas;
 *   betahilolo_d has space on the device for the second lowest doubles
 *             of the betas;
 *   betalololo_d has space on the device for the lowest doubles of the betas;
 *   add       current number of additions and subtractions;
 *   mul       current number of multiplications;
 *   div       current number of divisions;
 *   sqrtfun   current number of calls to sqrt();
 *   verbose   is the verbose flag.
 *
 * ON RETURN :
 *   vhihihi_h are the highest doubles of v,
 *             the next Householder vector on the host, if verbose;
 *   vlohihi_h are the second highest doubles of v;
 *   vhilohi_h are the third highest doubles of v;
 *   vlolohi_h are the fourth highest doubles of v;
 *   vhihilo_h are the fourth lowest doubles of v;
 *   vlohilo_h are the third lowest doubles of v;
 *   vhilolo_h are the second lowest doubles of v;
 *   vlololo_h are the lowest doubles of v;
 *   Vhihihi_d are the highest doubles of v, on the device
 *   Vlohihi_d are the second highest doubles of v, on the device;
 *   Vhilohi_d are the third highest doubles of v, on the device;
 *   Vlolohi_d are the fourth highest doubles of v, on the device;
 *   Vhihilo_d are the fourth lowest doubles of v, on the device;
 *   Vlohilo_d are the third lowest doubles of v, on the device;
 *   Vhilolo_d are the second lowest doubles of v, on the device;
 *   Vlololo_d are the lowest doubles of v, on the device;
 *   betahihihi_h has the updated vector of the highest doubles of the betas,
 *             if verbose;
 *   betalohihi_h has the updated vector of the second highest doubles
 *             of the betas, if verbose;
 *   betahilohi_h has the updated vector of the third highest doubles
 *             of the betas, if verbose;
 *   betalolohi_h has the updated vector of the fourth highest doubles
 *             of the betas, if verbose;
 *   betahihilo_h has the updated vector of the fourth lowest doubles
 *             of the betas, if verbose;
 *   betalohilo_h has the updated vector of the third lowest doubles
 *             of the betas, if verbose;
 *   betahilolo_h has the updated vector of the second lowest doubles
 *             of the betas, if verbose;
 *   betalololo_h has the updated vector of the lowest doubles of the betas,
 *             if verbose;
 *   betahihihi_d has the highest double of the next beta constant;
 *   betalohihi_d has the second highest double of the next beta constant;
 *   betahilohi_d has the third highest double of the next beta constant;
 *   betalolohi_d has the fourth highest double of the next beta constant;
 *   betahihilo_d has the fourth lowest double of the next beta constant;
 *   betalohilo_d has the third lowest double of the next beta constant;
 *   betahilolo_d has the second lowest double of the next beta constant;
 *   betalololo_d has the lowest double of the next beta constant;
 *   lapms     elapsed time spent by the kernel;
 *   add       accumulated number of additions and subtractions;
 *   mul       accumulated number of multiplications;
 *   div       accumulated number of divisions;
 *   sqrtfun   accumulated number of calls to sqrt(). */

void GPU_cmplx8_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *vrehihihi_h, double *vrelohihi_h,
   double *vrehilohi_h, double *vrelolohi_h,
   double *vrehihilo_h, double *vrelohilo_h,
   double *vrehilolo_h, double *vrelololo_h,
   double *vimhihihi_h, double *vimlohihi_h,
   double *vimhilohi_h, double *vimlolohi_h,
   double *vimhihilo_h, double *vimlohilo_h,
   double *vimhilolo_h, double *vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
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

void GPU_dbl8_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *vhihihi_h, double *vlohihi_h, double *vhilohi_h, double *vlolohi_h,
   double *vhihilo_h, double *vlohilo_h, double *vhilolo_h, double *vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *sumshihihi_h, double *sumslohihi_h,
   double *sumshilohi_h, double *sumslolohi_h,
   double *sumshihilo_h, double *sumslohilo_h,
   double *sumshilolo_h, double *sumslololo_h,
   double *sumshihihi_d, double *sumslohihi_d,
   double *sumshilohi_d, double *sumslolohi_d,
   double *sumshihilo_d, double *sumslohilo_d,
   double *sumshilolo_d, double *sumslololo_d,
   double *sigmahihihi_h, double *sigmalohihi_h, 
   double *sigmahilohi_h, double *sigmalolohi_h, 
   double *sigmahihilo_h, double *sigmalohilo_h, 
   double *sigmahilolo_h, double *sigmalololo_h, 
   double *sigmahihihi_d, double *sigmalohihi_d,
   double *sigmahilohi_d, double *sigmalolohi_d,
   double *sigmahihilo_d, double *sigmalohilo_d,
   double *sigmahilolo_d, double *sigmalololo_d,
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

void GPU_cmplx8_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *vrehihihi_h, double *vrelohihi_h,
   double *vrehilohi_h, double *vrelolohi_h,
   double *vrehihilo_h, double *vrelohilo_h,
   double *vrehilolo_h, double *vrelololo_h,
   double *vimhihihi_h, double *vimlohihi_h,
   double *vimhilohi_h, double *vimlolohi_h,
   double *vimhihilo_h, double *vimlohilo_h,
   double *vimhilolo_h, double *vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *sumshihihi_h, double *sumslohihi_h,
   double *sumshilohi_h, double *sumslolohi_h,
   double *sumshihilo_h, double *sumslohilo_h,
   double *sumshilolo_h, double *sumslololo_h,
   double *sumshihihi_d, double *sumslohihi_d,
   double *sumshilohi_d, double *sumslolohi_d,
   double *sumshihilo_d, double *sumslohilo_d,
   double *sumshilolo_d, double *sumslololo_d,
   double *sigmahihihi_h, double *sigmalohihi_h,
   double *sigmahilohi_h, double *sigmalolohi_h,
   double *sigmahihilo_h, double *sigmalohilo_h,
   double *sigmahilolo_h, double *sigmalololo_h,
   double *sigmahihihi_d, double *sigmalohihi_d,
   double *sigmahilohi_d, double *sigmalolohi_d,
   double *sigmahihilo_d, double *sigmalohilo_d,
   double *sigmahilolo_d, double *sigmalololo_d,
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

void GPU_dbl8_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
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

void GPU_cmplx8_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
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

void GPU_dbl8_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *RTdotvhihihi_h, double *RTdotvlohihi_h,
   double *RTdotvhilohi_h, double *RTdotvlolohi_h,
   double *RTdotvhihilo_h, double *RTdotvlohilo_h,
   double *RTdotvhilolo_h, double *RTdotvlololo_h,
   double *RTdotvhihihi_d, double *RTdotvlohihi_d,
   double *RTdotvhilohi_d, double *RTdotvlolohi_d,
   double *RTdotvhihilo_d, double *RTdotvlohilo_d,
   double *RTdotvhilolo_d, double *RTdotvlololo_d,
   double *whihihi_h, double *wlohihi_h, double *whilohi_h, double *wlolohi_h,
   double *whihilo_h, double *wlohilo_h, double *whilolo_h, double *wlololo_h,
   double *whihihi_d, double *wlohihi_d, double *whilohi_d, double *wlolohi_d,
   double *whihilo_d, double *wlohilo_d, double *whilolo_d, double *wlololo_d,
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

void GPU_cmplx8_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *RHdotvrehihihi_h, double *RHdotvrelohihi_h,
   double *RHdotvrehilohi_h, double *RHdotvrelolohi_h,
   double *RHdotvrehihilo_h, double *RHdotvrelohilo_h,
   double *RHdotvrehilolo_h, double *RHdotvrelololo_h,
   double *RHdotvimhihihi_h, double *RHdotvimlohihi_h,
   double *RHdotvimhilohi_h, double *RHdotvimlolohi_h,
   double *RHdotvimhihilo_h, double *RHdotvimlohilo_h,
   double *RHdotvimhilolo_h, double *RHdotvimlololo_h,
   double *RHdotvrehihihi_d, double *RHdotvrelohihi_d,
   double *RHdotvrehilohi_d, double *RHdotvrelolohi_d,
   double *RHdotvrehihilo_d, double *RHdotvrelohilo_d,
   double *RHdotvrehilolo_d, double *RHdotvrelololo_d,
   double *RHdotvimhihihi_d, double *RHdotvimlohihi_d,
   double *RHdotvimhilohi_d, double *RHdotvimlolohi_d,
   double *RHdotvimhihilo_d, double *RHdotvimlohilo_d,
   double *RHdotvimhilolo_d, double *RHdotvimlololo_d,
   double *wrehihihi_h, double *wrelohihi_h,
   double *wrehilohi_h, double *wrelolohi_h,
   double *wrehihilo_h, double *wrelohilo_h,
   double *wrehilolo_h, double *wrelololo_h,
   double *wimhihihi_h, double *wimlohihi_h,
   double *wimhilohi_h, double *wimlolohi_h,
   double *wimhihilo_h, double *wimlohilo_h,
   double *wimhilolo_h, double *wimlololo_h,
   double *wrehihihi_d, double *wrelohihi_d,
   double *wrehilohi_d, double *wrelolohi_d,
   double *wrehihilo_d, double *wrelohilo_d,
   double *wrehilolo_d, double *wrelololo_d,
   double *wimhihihi_d, double *wimlohihi_d,
   double *wimhilohi_d, double *wimlolohi_d,
   double *wimhihilo_d, double *wimlohilo_d,
   double *wimhilolo_d, double *wimlololo_d,
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

void GPU_dbl8_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vhihihi_h, double *Vlohihi_h, double *Vhilohi_h, double *Vlolohi_h,
   double *Vhihilo_h, double *Vlohilo_h, double *Vhilolo_h, double *Vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *Whihihi_h, double *Wlohihi_h, double *Whilohi_h, double *Wlolohi_h,
   double *Whihilo_h, double *Wlohilo_h, double *Whilolo_h, double *Wlololo_h,
   double *Whihihi_d, double *Wlohihi_d, double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d, double *Whilolo_d, double *Wlololo_d,
   double *WYThihihi_h, double *WYTlohihi_h,
   double *WYThilohi_h, double *WYTlolohi_h,
   double *WYThihilo_h, double *WYTlohilo_h,
   double *WYThilolo_h, double *WYTlololo_h,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
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

void GPU_cmplx8_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vrehihihi_h, double *Vrelohihi_h,
   double *Vrehilohi_h, double *Vrelolohi_h,
   double *Vrehihilo_h, double *Vrelohilo_h,
   double *Vrehilolo_h, double *Vrelololo_h,
   double *Vimhihihi_h, double *Vimlohihi_h,
   double *Vimhilohi_h, double *Vimlolohi_h,
   double *Vimhihilo_h, double *Vimlohilo_h,
   double *Vimhilolo_h, double *Vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d, 
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *Wrehihihi_h, double *Wrelohihi_h,
   double *Wrehilohi_h, double *Wrelolohi_h,
   double *Wrehihilo_h, double *Wrelohilo_h,
   double *Wrehilolo_h, double *Wrelololo_h,
   double *Wimhihihi_h, double *Wimlohihi_h,
   double *Wimhilohi_h, double *Wimlolohi_h,
   double *Wimhihilo_h, double *Wimlohilo_h,
   double *Wimhilolo_h, double *Wimlololo_h,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *WYHrehihihi_h, double *WYHrelohihi_h,
   double *WYHrehilohi_h, double *WYHrelolohi_h,
   double *WYHrehihilo_h, double *WYHrelohilo_h,
   double *WYHrehilolo_h, double *WYHrelololo_h,
   double *WYHimhihihi_h, double *WYHimlohihi_h,
   double *WYHimhilohi_h, double *WYHimlolohi_h,
   double *WYHimhihilo_h, double *WYHimlohilo_h,
   double *WYHimhilolo_h, double *WYHimlololo_h,
   double *WYHrehihihi_d, double *WYHrelohihi_d,
   double *WYHrehilohi_d, double *WYHrelolohi_d,
   double *WYHrehihilo_d, double *WYHrelohilo_d,
   double *WYHrehilolo_d, double *WYHrelololo_d,
   double *WYHimhihihi_d, double *WYHimlohihi_d,
   double *WYHimhilohi_d, double *WYHimlolohi_d,
   double *WYHimhihilo_d, double *WYHimlohilo_d,
   double *WYHimhilolo_d, double *WYHimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the W in the WY representation,
 *   on complex data.
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
 *   Vrehi_h  high doubles of the real parts
 *            of the Householder vectors, if verbose;
 *   Vrelo_h  low doubles of the real parts
 *            of the Householder vectors, if verbose;
 *   Vimhi_h  high doubles of the imaginary parts
 *            of the Householder vectors, if verbose;
 *   Vimlo_h  low doubles of the imaginary parts
 *            of the Householder vectors, if verbose;
 *   Vrehi_d  high doubles of the the real parts of V, on the device;
 *   Vrelo_d  low doubles of the real parts of V, on the device;
 *   Vimhi_d  high doubles of the imaginary parts of V, on the device;
 *   Vimlo_d  low doubles of the imaginary parts of V, on the device;
 *   Wrehi_h  has space for the high doubles of the real parts of W,
 *            if verbose;
 *   Wrelo_h  has space for the low doubles of the real parts of W,
 *            if verbose;
 *   Wimhi_h  has space for the high doubles of the imaginary parts of W,
 *            if verbose;
 *   Wimlo_h  has space for the low doubles of the imaginary parts of W,
 *            if verbose;
 *   Wrehi_d  has space for the high doubles of the real parts of W,
 *            on the device;
 *   Wrelo_d  has space for the low doubles of the real parts of W,
 *            on the device;
 *   Wimhi_d  has space for the high doubles of the imaginary parts of W,
 *            on the device;
 *   Wimlo_d  has space for the low doubles of the imaginary parts of W,
 *            on the device;
 *   WYHrehi_h has space for the high doubles of the real parts of W*Y^H,
 *            if verbose;
 *   WYHrelo_h has space for the low doubles of the real parts of W*Y^H,
 *            if verbose;
 *   WYHimhi_h has space for the high doubles of the imaginary parts of W*Y^H,
 *            if verbose;
 *   WYHimlo_h has space for the low doubles of the imaginary parts of W*Y^H,
 *            if verbose;
 *   WYHrehi_d has space for the high doubles of the real parts of W*Y^H,
 *            on the device;
 *   WYHrelo_d has space for the low doubles of the real parts of W*Y^H,
 *            on the device;
 *   WYHimhi_d has space for the high doubles of the imaginary parts of W*Y^H,
 *            on the device;
 *   WYHimlo_d has space for the low doubles of the imaginary parts of W*Y^H,
 *            on the device;
 *   betahi_h has space for the betas if verbose;
 *   betalo_h has space for the betas if verbose;
 *   betahi_d has space on the device for the betas;
 *   betalo_d has space on the device for the betas;
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

void GPU_dbl8_small_WYT
 ( int nrows, int szt,
   double *Whihihi_d, double *Wlohihi_d,
   double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d,
   double *Whilolo_d, double *Wlololo_d,
   double *Yhihihi_d, double *Ylohihi_d,
   double *Yhilohi_d, double *Ylolohi_d,
   double *Yhihilo_d, double *Ylohilo_d,
   double *Yhilolo_d, double *Ylololo_d,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *WYThihihi_h, double *WYTlohihi_h,
   double *WYThilohi_h, double *WYTlolohi_h,
   double *WYThihilo_h, double *WYTlohilo_h,
   double *WYThilolo_h, double *WYTlololo_h,
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

void GPU_cmplx8_small_WYH
 ( int nrows, int szt,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *Yrehihihi_d, double *Yrelohihi_d,
   double *Yrehilohi_d, double *Yrelolohi_d,
   double *Yrehihilo_d, double *Yrelohilo_d,
   double *Yrehilolo_d, double *Yrelololo_d,
   double *Yimhihihi_d, double *Yimlohihi_d,
   double *Yimhilohi_d, double *Yimlolohi_d,
   double *Yimhihilo_d, double *Yimlohilo_d,
   double *Yimhilolo_d, double *Yimlololo_d,
   double *WYTrehihihi_d, double *WYTrelohihi_d,
   double *WYTrehilohi_d, double *WYTrelolohi_d,
   double *WYTrehihilo_d, double *WYTrelohilo_d,
   double *WYTrehilolo_d, double *WYTrelololo_d,
   double *WYTimhihihi_d, double *WYTimlohihi_d,
   double *WYTimhilohi_d, double *WYTimlolohi_d,
   double *WYTimhihilo_d, double *WYTimlohilo_d,
   double *WYTimhilolo_d, double *WYTimlololo_d,
   double *WYTrehihihi_h, double *WYTrelohihi_h,
   double *WYTrehilohi_h, double *WYTrelolohi_h,
   double *WYTrehihilo_h, double *WYTrelohilo_h,
   double *WYTrehilolo_h, double *WYTrelololo_h,
   double *WYTimhihihi_h, double *WYTimlohihi_h,
   double *WYTimhilohi_h, double *WYTimlolohi_h,
   double *WYTimhihilo_h, double *WYTimlohilo_h,
   double *WYTimhilolo_h, double *WYTimlololo_h,
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

void GPU_dbl8_small_YWT
 ( int nrows, int szt, int idx,
   double *Yhihihi_d, double *Ylohihi_d, double *Yhilohi_d, double *Ylolohi_d,
   double *Yhihilo_d, double *Ylohilo_d, double *Yhilolo_d, double *Ylololo_d,
   double *Whihihi_d, double *Wlohihi_d, double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d, double *Whilolo_d, double *Wlololo_d,
   double *YWThihihi_d, double *YWTlohihi_d,
   double *YWThilohi_d, double *YWTlolohi_d,
   double *YWThihilo_d, double *YWTlohilo_d,
   double *YWThilolo_d, double *YWTlololo_d,
   double *YWThihihi_h, double *YWTlohihi_h,
   double *YWThilohi_h, double *YWTlolohi_h,
   double *YWThihilo_h, double *YWTlohilo_h,
   double *YWThilolo_h, double *YWTlololo_h,
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

void GPU_cmplx8_small_YWH
 ( int nrows, int szt, int idx,
   double *Yrehihihi_d, double *Yrelohihi_d,
   double *Yrehilohi_d, double *Yrelolohi_d,
   double *Yrehihilo_d, double *Yrelohilo_d,
   double *Yrehilolo_d, double *Yrelololo_d,
   double *Yimhihihi_d, double *Yimlohihi_d,
   double *Yimhilohi_d, double *Yimlolohi_d,
   double *Yimhihilo_d, double *Yimlohilo_d,
   double *Yimhilolo_d, double *Yimlololo_d,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *YWTrehihihi_d, double *YWTrelohihi_d,
   double *YWTrehilohi_d, double *YWTrelolohi_d,
   double *YWTrehihilo_d, double *YWTrelohilo_d,
   double *YWTrehilolo_d, double *YWTrelololo_d,
   double *YWTimhihihi_d, double *YWTimlohihi_d,
   double *YWTimhilohi_d, double *YWTimlolohi_d,
   double *YWTimhihilo_d, double *YWTimlohilo_d,
   double *YWTimhilolo_d, double *YWTimlololo_d,
   double *YWTrehihihi_h, double *YWTrelohihi_h,
   double *YWTrehilohi_h, double *YWTrelolohi_h,
   double *YWTrehihilo_h, double *YWTrelohilo_h,
   double *YWTrehilolo_h, double *YWTrelololo_h,
   double *YWTimhihihi_h, double *YWTimlohihi_h,
   double *YWTimhilohi_h, double *YWTimlolohi_h,
   double *YWTimhihilo_h, double *YWTimlohilo_h,
   double *YWTimhilolo_h, double *YWTimlololo_h,
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

void GPU_dbl8_small_QWYT
 ( int dim, int szt, int idx,
   double *Qhihihi_d, double *Qlohihi_d, double *Qhilohi_d, double *Qlolohi_d,
   double *Qhihilo_d, double *Qlohilo_d, double *Qhilolo_d, double *Qlololo_d,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *QWYThihihi_d, double *QWYTlohihi_d,
   double *QWYThilohi_d, double *QWYTlolohi_d,
   double *QWYThihilo_d, double *QWYTlohilo_d,
   double *QWYThilolo_d, double *QWYTlololo_d,
   double *QWYThihihi_h, double *QWYTlohihi_h,
   double *QWYThilohi_h, double *QWYTlolohi_h,
   double *QWYThihilo_h, double *QWYTlohilo_h,
   double *QWYThilolo_h, double *QWYTlololo_h,
   double *Qhihihi_h, double *Qlohihi_h, double *Qhilohi_h, double *Qlolohi_h,
   double *Qhihilo_h, double *Qlohilo_h, double *Qhilolo_h, double *Qlololo_h,
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

void GPU_cmplx8_small_QWYH
 ( int dim, int szt, int idx,
   double *Qrehihihi_d, double *Qrelohihi_d,
   double *Qrehilohi_d, double *Qrelolohi_d,
   double *Qrehihilo_d, double *Qrelohilo_d,
   double *Qrehilolo_d, double *Qrelololo_d,
   double *Qimhihihi_d, double *Qimlohihi_d,
   double *Qimhilohi_d, double *Qimlolohi_d,
   double *Qimhihilo_d, double *Qimlohilo_d,
   double *Qimhilolo_d, double *Qimlololo_d,
   double *WYTrehihihi_d, double *WYTrelohihi_d,
   double *WYTrehilohi_d, double *WYTrelolohi_d,
   double *WYTrehihilo_d, double *WYTrelohilo_d,
   double *WYTrehilolo_d, double *WYTrelololo_d,
   double *WYTimhihihi_d, double *WYTimlohihi_d,
   double *WYTimhilohi_d, double *WYTimlolohi_d,
   double *WYTimhihilo_d, double *WYTimlohilo_d,
   double *WYTimhilolo_d, double *WYTimlololo_d,
   double *QWYTrehihihi_d, double *QWYTrelohihi_d,
   double *QWYTrehilohi_d, double *QWYTrelolohi_d,
   double *QWYTrehihilo_d, double *QWYTrelohilo_d,
   double *QWYTrehilolo_d, double *QWYTrelololo_d,
   double *QWYTimhihihi_d, double *QWYTimlohihi_d,
   double *QWYTimhilohi_d, double *QWYTimlolohi_d,
   double *QWYTimhihilo_d, double *QWYTimlohilo_d,
   double *QWYTimhilolo_d, double *QWYTimlololo_d,
   double *QWYTrehihihi_h, double *QWYTrelohihi_h,
   double *QWYTrehilohi_h, double *QWYTrelolohi_h,
   double *QWYTrehihilo_h, double *QWYTrelohilo_h,
   double *QWYTrehilolo_h, double *QWYTrelololo_h,
   double *QWYTimhihihi_h, double *QWYTimlohihi_h,
   double *QWYTimhilohi_h, double *QWYTimlolohi_h,
   double *QWYTimhihilo_h, double *QWYTimlohilo_h,
   double *QWYTimhilolo_h, double *QWYTimlololo_h,
   double *Qrehihihi_h, double *Qrelohihi_h,
   double *Qrehilohi_h, double *Qrelolohi_h,
   double *Qrehihilo_h, double *Qrelohilo_h,
   double *Qrehilolo_h, double *Qrelololo_h,
   double *Qimhihihi_h, double *Qimlohihi_h,
   double *Qimhilohi_h, double *Qimlolohi_h,
   double *Qimhihilo_h, double *Qimlohilo_h,
   double *Qimhilolo_h, double *Qimlololo_h,
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

void GPU_dbl8_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWThihihi_d, double *YWTlohihi_d,
   double *YWThilohi_d, double *YWTlolohi_d,
   double *YWThihilo_d, double *YWTlohilo_d,
   double *YWThilolo_d, double *YWTlololo_d,
   double *Chihihi_d, double *Clohihi_d, double *Chilohi_d, double *Clolohi_d,
   double *Chihilo_d, double *Clohilo_d, double *Chilolo_d, double *Clololo_d,
   double *YWTChihihi_d, double *YWTClohihi_d,
   double *YWTChilohi_d, double *YWTClolohi_d,
   double *YWTChihilo_d, double *YWTClohilo_d,
   double *YWTChilolo_d, double *YWTClololo_d,
   double *YWTChihihi_h, double *YWTClohihi_h,
   double *YWTChilohi_h, double *YWTClolohi_h,
   double *YWTChihilo_h, double *YWTClohilo_h,
   double *YWTChilolo_h, double *YWTClololo_h,
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

void GPU_cmplx8_small_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTrehihihi_d, double *YWTrelohihi_d,
   double *YWTrehilohi_d, double *YWTrelolohi_d,
   double *YWTrehihilo_d, double *YWTrelohilo_d,
   double *YWTrehilolo_d, double *YWTrelololo_d,
   double *YWTimhihihi_d, double *YWTimlohihi_d,
   double *YWTimhilohi_d, double *YWTimlolohi_d,
   double *YWTimhihilo_d, double *YWTimlohilo_d,
   double *YWTimhilolo_d, double *YWTimlololo_d,
   double *Crehihihi_d, double *Crelohihi_d,
   double *Crehilohi_d, double *Crelolohi_d,
   double *Crehihilo_d, double *Crelohilo_d,
   double *Crehilolo_d, double *Crelololo_d,
   double *Cimhihihi_d, double *Cimlohihi_d,
   double *Cimhilohi_d, double *Cimlolohi_d,
   double *Cimhihilo_d, double *Cimlohilo_d,
   double *Cimhilolo_d, double *Cimlololo_d,
   double *YWTCrehihihi_d, double *YWTCrelohihi_d,
   double *YWTCrehilohi_d, double *YWTCrelolohi_d,
   double *YWTCrehihilo_d, double *YWTCrelohilo_d,
   double *YWTCrehilolo_d, double *YWTCrelololo_d,
   double *YWTCimhihihi_d, double *YWTCimlohihi_d,
   double *YWTCimhilohi_d, double *YWTCimlolohi_d,
   double *YWTCimhihilo_d, double *YWTCimlohilo_d,
   double *YWTCimhilolo_d, double *YWTCimlololo_d,
   double *YWTCrehihihi_h, double *YWTCrelohihi_h,
   double *YWTCrehilohi_h, double *YWTCrelolohi_h,
   double *YWTCrehihilo_h, double *YWTCrelohilo_h,
   double *YWTCrehilolo_h, double *YWTCrelololo_h,
   double *YWTCimhihihi_h, double *YWTCimlohihi_h,
   double *YWTCimhilohi_h, double *YWTCimlolohi_h,
   double *YWTCimhihilo_h, double *YWTCimlohilo_h,
   double *YWTCimhilolo_h, double *YWTCimlololo_h,
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

void GPU_dbl8_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qhihihi_d, double *Qlohihi_d, double *Qhilohi_d, double *Qlolohi_d,
   double *Qhihilo_d, double *Qlohilo_d, double *Qhilolo_d, double *Qlololo_d,
   double *QWYThihihi_d, double *QWYTlohihi_d,
   double *QWYThilohi_d, double *QWYTlolohi_d,
   double *QWYThihilo_d, double *QWYTlohilo_d,
   double *QWYThilolo_d, double *QWYTlololo_d,
   double *Qhihihi_h, double *Qlohihi_h, double *Qhilohi_h, double *Qlolohi_h,
   double *Qhihilo_h, double *Qlohilo_h, double *Qhilolo_h, double *Qlololo_h,
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

void GPU_cmplx8_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qrehihihi_d, double *Qrelohihi_d,
   double *Qrehilohi_d, double *Qrelolohi_d,
   double *Qrehihilo_d, double *Qrelohilo_d,
   double *Qrehilolo_d, double *Qrelololo_d,
   double *Qimhihihi_d, double *Qimlohihi_d,
   double *Qimhilohi_d, double *Qimlolohi_d,
   double *Qimhihilo_d, double *Qimlohilo_d,
   double *Qimhilolo_d, double *Qimlololo_d,
   double *QWYTrehihihi_d, double *QWYTrelohihi_d,
   double *QWYTrehilohi_d, double *QWYTrelolohi_d,
   double *QWYTrehihilo_d, double *QWYTrelohilo_d,
   double *QWYTrehilolo_d, double *QWYTrelololo_d,
   double *QWYTimhihihi_d, double *QWYTimlohihi_d,
   double *QWYTimhilohi_d, double *QWYTimlolohi_d,
   double *QWYTimhihilo_d, double *QWYTimlohilo_d,
   double *QWYTimhilolo_d, double *QWYTimlololo_d,
   double *Qrehihihi_h, double *Qrelohihi_h,
   double *Qrehilohi_h, double *Qrelolohi_h,
   double *Qrehihilo_h, double *Qrelohilo_h,
   double *Qrehilolo_h, double *Qrelololo_h,
   double *Qimhihihi_h, double *Qimlohihi_h,
   double *Qimhilohi_h, double *Qimlolohi_h,
   double *Qimhihilo_h, double *Qimlohilo_h,
   double *Qimhilolo_h, double *Qimlololo_h,
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

void GPU_dbl8_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rhihihi_d, double *Rlohihi_d, double *Rhilohi_d, double *Rlolohi_d,
   double *Rhihilo_d, double *Rlohilo_d, double *Rhilolo_d, double *Rlololo_d,
   double *YWTChihihi_d, double *YWTClohihi_d,
   double *YWTChilohi_d, double *YWTClolohi_d,
   double *YWTChihilo_d, double *YWTClohilo_d,
   double *YWTChilolo_d, double *YWTClololo_d,
   double *Rhihihi_h, double *Rlohihi_h, double *Rhilohi_h, double *Rlolohi_h,
   double *Rhihilo_h, double *Rlohilo_h, double *Rhilolo_h, double *Rlololo_h,
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

void GPU_cmplx8_small_R_add_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *Rrehihihi_d, double *Rrelohihi_d,
   double *Rrehilohi_d, double *Rrelolohi_d,
   double *Rrehihilo_d, double *Rrelohilo_d,
   double *Rrehilolo_d, double *Rrelololo_d,
   double *Rimhihihi_d, double *Rimlohihi_d,
   double *Rimhilohi_d, double *Rimlolohi_d,
   double *Rimhihilo_d, double *Rimlohilo_d,
   double *Rimhilolo_d, double *Rimlololo_d,
   double *YWTCrehihihi_d, double *YWTCrelohihi_d,
   double *YWTCrehilohi_d, double *YWTCrelolohi_d,
   double *YWTCrehihilo_d, double *YWTCrelohilo_d,
   double *YWTCrehilolo_d, double *YWTCrelololo_d,
   double *YWTCimhihihi_d, double *YWTCimlohihi_d,
   double *YWTCimhilohi_d, double *YWTCimlolohi_d,
   double *YWTCimhihilo_d, double *YWTCimlohilo_d,
   double *YWTCimhilolo_d, double *YWTCimlololo_d,
   double *Rrehihihi_h, double *Rrelohihi_h,
   double *Rrehilohi_h, double *Rrelolohi_h,
   double *Rrehihilo_h, double *Rrelohilo_h,
   double *Rrehilolo_h, double *Rrelololo_h,
   double *Rimhihihi_h, double *Rimlohihi_h,
   double *Rimhilohi_h, double *Rimlolohi_h,
   double *Rimhihilo_h, double *Rimlohilo_h,
   double *Rimhilolo_h, double *Rimlololo_h,
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

void GPU_dbl8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
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
 *   Ahihihi  highest doubles of an nrows-by-ncols matrix A,
 *            stored as nrows arrays of ncols numbers;
 *   Alohihi  second highest doubles of A;
 *   Ahilohi  third highest doubles of A;
 *   Alolohi  fourth highest doubles of A;
 *   Ahihilo  fourth lowest doubles of A;
 *   Alohilo  third lowest doubles of A;
 *   Ahilolo  second lowest doubles of A;
 *   Alololo  lowest doubles of A;
 *   Qhihihi  space for an nrows-by-nrows matrix;
 *   Qlohihi  space for an nrows-by-nrows matrix;
 *   Qhilohi  space for an nrows-by-nrows matrix;
 *   Qlolohi  space for an nrows-by-nrows matrix;
 *   Qhihilo  space for an nrows-by-nrows matrix;
 *   Qlohilo  space for an nrows-by-nrows matrix;
 *   Qhilolo  space for an nrows-by-nrows matrix;
 *   Qlololo  space for an nrows-by-nrows matrix;
 *   Rhihihi  space for an nrows-by-ncols matrix;
 *   Rlohihi  space for an nrows-by-ncols matrix;
 *   Rhilohi  space for an nrows-by-ncols matrix;
 *   Rlolohi  space for an nrows-by-ncols matrix;
 *   Rhihilo  space for an nrows-by-ncols matrix;
 *   Rlohilo  space for an nrows-by-ncols matrix;
 *   Rhilolo  space for an nrows-by-ncols matrix;
 *   Rlololo  space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhihihi  highest doubles of an orthogonal matrix,
              transpose(Q)*A = R;
 *   Qlohihi  second highest doubles of Q;
 *   Qhilohi  third highest doubles of Q;
 *   Qlolohi  fourth highest doubles of Q;
 *   Qhihilo  fourth lowest doubles of Q;
 *   Qlohilo  third lowest doubles of Q;
 *   Qhilolo  second lowest doubles of Q;
 *   Qlololo  lowest doubles of Q;
 *   Rhihihi  highest doubles of the reduced upper triangular form R of A;
 *   Rlohihi  second highest doubles of R;
 *   Rhilohi  third highest doubles of R;
 *   Rlolohi  fourth highest doubles of R;
 *   Rhihilo  fourth lowest doubles of R;
 *   Rlohilo  third lowest doubles of R;
 *   Rhilolo  second lowest doubles of R;
 *   Rlololo  lowest doubles of R;
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

void GPU_cmplx8_blocked_houseqr
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
 *   nrows     number of rows of A;
 *   ncols     number of columns of A;
 *   szt       size of each block;
 *   nbt       number of tiles, ncols = szt*nbt;
 *   Arehihihi are the highest doubles of the real parts of A,
 *             stored as nrows arrays of ncols numbers;
 *   Arelohihi are the second highest doubles of the real parts of A;
 *   Arehilohi are the third highest doubles of the real parts of A;
 *   Arelolohi are the fourth highest doubles of the real parts of A;
 *   Arehihilo are the fourth lowest doubles of the real parts of A;
 *   Arehihilo are the third lowest doubles of the real parts of A;
 *   Arelohilo are the second lowest doubles of the real parts of A;
 *   Arelololo are the lowest doubles of the real parts of A;
 *   Aimhihihi are the highest doubles of the imaginary parts of A;
 *   Aimlohihi are the second highest doubles of the imaginary parts of A;
 *   Aimhilohi are the third highest doubles of the imaginary parts of A;
 *   Aimlolohi are the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo are the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo are the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo are the second lowest doubles of the imaginary parts of A;
 *   Aimlololo are the lowest doubles of the imaginary parts of A;
 *   Qrehihihi has space for an nrows-by-nrows matrix;
 *   Qrelohihi has space for an nrows-by-nrows matrix;
 *   Qrehilohi has space for an nrows-by-nrows matrix;
 *   Qrelolohi has space for an nrows-by-nrows matrix;
 *   Qrehihilo has space for an nrows-by-nrows matrix;
 *   Qrelohilo has space for an nrows-by-nrows matrix;
 *   Qrehilolo has space for an nrows-by-nrows matrix;
 *   Qrelololo has space for an nrows-by-nrows matrix;
 *   Qimhihihi has space for an nrows-by-nrows matrix;
 *   Qimlohihi has space for an nrows-by-nrows matrix;
 *   Qimhilohi has space for an nrows-by-nrows matrix;
 *   Qimlolohi has space for an nrows-by-nrows matrix;
 *   Qimhihilo has space for an nrows-by-nrows matrix;
 *   Qimlohilo has space for an nrows-by-nrows matrix;
 *   Qimhilolo has space for an nrows-by-nrows matrix;
 *   Qimlololo has space for an nrows-by-nrows matrix;
 *   Rrehihihi has space for an nrows-by-ncols matrix;
 *   Rrelohihi has space for an nrows-by-ncols matrix;
 *   Rrehilohi has space for an nrows-by-ncols matrix;
 *   Rrelolohi has space for an nrows-by-ncols matrix;
 *   Rrehihilo has space for an nrows-by-ncols matrix;
 *   Rrelohilo has space for an nrows-by-ncols matrix;
 *   Rrehilolo has space for an nrows-by-ncols matrix;
 *   Rrelololo has space for an nrows-by-ncols matrix;
 *   Rimhihihi has space for an nrows-by-ncols matrix;
 *   Rimlohihi has space for an nrows-by-ncols matrix;
 *   Rimhilohi has space for an nrows-by-ncols matrix;
 *   Rimlolohi has space for an nrows-by-ncols matrix;
 *   Rimhihilo has space for an nrows-by-ncols matrix;
 *   Rimlohilo has space for an nrows-by-ncols matrix;
 *   Rimhilolo has space for an nrows-by-ncols matrix;
 *   Rimlololo has space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qrehihihi has the highest doubles of the real parts
 *             of an orthogonal matrix Q, transpose(Q)*A = R;
 *   Qrelohihi has the second highest doubles of the real parts of Q;
 *   Qrehilohi has the third highest doubles of the real parts of Q;
 *   Qrelolohi has the fourth highest doubles of the real parts of Q;
 *   Qrehihilo has the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo has the third lowest doubles of the real parts of Q;
 *   Qrehilolo has the second lowest doubles of the real parts of Q;
 *   Qrelololo has the lowest doubles of the real parts of Q;
 *   Qimhihihi has the highest doubles of the imaginary parts of Q;
 *   Qimlohihi has the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi has the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi has the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo has the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo has the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo has the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo has the lowest doubles of the imaginary parts of Q;
 *   Rrehihihi has the highest doubles of the real parts of the reduced
 *             upper triangular form R of A;
 *   Rrelohihi has the second highest doubles of the real parts of R;
 *   Rrehilohi has the third highest doubles of the real parts of R;
 *   Rrelolohi has the fourth highest doubles of the real parts of R;
 *   Rrehihilo has the fourth lowest doubles of the real parts of R;
 *   Rrelohilo has the third lowest doubles of the real parts of R;
 *   Rrehilolo has the second lowest doubles of the real parts of R;
 *   Rrelololo has the lowest doubles of the real parts of R;
 *   Rimhihihi has the highest doubles of the imaginary parts of R;
 *   Rimlohihi has the second highest doubles of the imaginary parts of R;
 *   Rimhilohi has the third highest doubles of the imaginary parts of R;
 *   Rimlolohi has the fourth highest doubles of the imaginary parts of R;
 *   Rimhihilo has the fourth lowest doubles of the imaginary parts of R;
 *   Rimlohilo has the third lowest doubles of the imaginary parts of R;
 *   Rimhilolo has the second lowest doubles of the imaginary parts of R;
 *   Rimlololo has the lowest doubles of the imaginary parts of R;
 *   houselapms is the elapsed time spent by the kernel
 *             to compute the Householder vector and the beta;
 *   RHvlapms  is the elapsed time spent to compute beta*R^H*v;
 *   tileRlapms is the elapsed time spent by the kernel
 *             to reduce one tile;
 *   vb2Wlapms is the elapsed time spent by the kernel
 *             to compute the W representation;
 *   WYHlapms  is the elapsed time spent by the kernel
 *             to compute the W*Y^H product;
 *   QWYHlapms is the elapsed time spent by the kernel
 *             to compute the Q*WYH product;
 *   Qaddlapms is the elapsed time spent by the kernel
 *             to compute Q by adding the Q*W*Y^H matrix;
 *   YWHlapms  is the elapsed time spent by the kernel
 *             to compute the Y*W^H product;
 *   YWHClapms is the elapsed time spent by the kernel
 *             to compute the YWH*C product;
 *   Raddlapms is the elapsed time spent by the kernel
 *             to compute R by adding the Y*W^H*C matrix;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt    counts the number of additions and subtractions;
 *   mulcnt    counts the number of multiplications;
 *   divcnt    counts the number of divisions;
 *   sqrtcnt   counts the number of calls to sqrt(). */

#endif
