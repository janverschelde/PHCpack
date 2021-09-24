/* The file dbl_baqr_flopcounts.cpp defines the functions with prototypes in
 * the file dbl_baqr_flopcounts.h. */

#include "dbl_baqr_flopcounts.h"

void flopcount_dbl_small_house
 ( int dim, int dimLog2, long long int *add, long long int *mul,
   long long int *div, long long int *sqrtfun )
{
   *add += dimLog2 + 4;
   *mul += dim + 3;
   *div += dim + 1;
   *sqrtfun += 1;
}

void flopcount_cmplx_small_house
 ( int dim, int dimLog2, long long int *add,
   long long int *mul, long long int *div, long long int *sqrtfun )
{
   *add += 3*dim + dimLog2 + 6;
   *mul += 6*dim + 6;
   *div += 2*dim + 2;
   *sqrtfun += 3;
}

void flopcount_dbl_large_sum_of_squares
 ( int nblocks, int szt, int sztLog2, long long int *add, long long int *mul )
{
   *add += nblocks*szt*sztLog2;
   *mul += nblocks*szt;
}

void flopcount_cmplx_large_sum_of_squares
 ( int nblocks, int szt, int sztLog2, long long int *add, long long int *mul )
{
   *add += nblocks*szt*(1+sztLog2);
   *mul += 2*nblocks*szt;
}

void flopcount_dbl_sum_accumulator
 ( int nbt, int nbtLog2, long long int *add )
{
   *add += nbt*nbtLog2;
}

void flopcount_dbl_normalize ( int nblocks, int szt, long long int *div )
{
   *div += nblocks*szt;
}

void flopcount_cmplx_normalize
 ( int nblocks, int szt, long long int *add, long long int *mul )
{
   *add += 2*nblocks*szt;
   *mul += 4*nblocks*szt;
}

void flopcount_dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul )
{
   const int nbthreads = ncols - k;

   *add += 2*(nrows - k)*nbthreads;
   *mul += 2*(nrows - k)*nbthreads + nbthreads;
}

void flopcount_cmplx_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul )
{
   const int nbthreads = ncols - k;

   *add += 8*(nrows - k)*nbthreads;
   *mul += 8*(nrows - k)*nbthreads + nbthreads;
}

void flopcount_dbl_small_betaRTv 
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul )
{
   const int nbthreads = nrows - k;

   *add += (nrows-k)*nbthreads;
   *mul += (nrows-k)*nbthreads + nbthreads;
}

void flopcount_cmplx_small_betaRHv 
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul )
{
   const int nbthreads = nrows - k;

   *add += 4*(nrows-k)*nbthreads;
   *mul += 4*(nrows-k)*nbthreads + 2*nbthreads;
}

void flopcount_dbl_RTdotv ( int nrows, int szt, long long int *mul )
{
   int nbthreads = nrows*szt;

   *mul += nbthreads;
}

void flopcount_cmplx_RHdotv
 ( int nrows, int szt, long long int *add, long long int *mul )
{
   int nbthreads = nrows*szt;

   *add += 2*nbthreads;
   *mul += 4*nbthreads;
}

void flopcount_dbl_sum_betaRTdotv
 ( int nrows, int dim, long long int *add, long long int *mul )
{
   *add += dim*nrows;
   *mul += dim;
}

void flopcount_cmplx_sum_betaRHdotv
 ( int nrows, int dim, long long int *add, long long int *mul )
{
   *add += 2*dim*nrows;
   *mul += 2*dim;
}

void flopcount_dbl_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul )
{
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   // total number of entries in R that will be modified
   // equals the number of threads that compute
   const int nbthreads = (nrows - k)*(endcol - k);

   *add += nbthreads;
   *mul += nbthreads;
}

void flopcount_cmplx_medium_subvbetaRHv
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul )
{
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   // total number of entries in R that will be modified
   // equals the number of threads that compute
   const int nbthreads = (nrows - k)*(endcol - k);

   *add += 4*nbthreads;
   *mul += 4*nbthreads;
}

void flopcount_dbl_VB_to_W
 ( int nrows, int ncols, long long int *add, long long int *mul )
{
   // the number of threads equals nrows
   *add += nrows*((ncols-1)*ncols/2*(1+nrows) + ncols-1);
   *mul += nrows*(1 + (ncols-1)*ncols + ncols-1);
}

void flopcount_cmplx_VB_to_W
 ( int nrows, int ncols, long long int *add, long long int *mul )
{
   // the number of threads equals nrows
   *add += nrows*((ncols-1)*ncols/2*(6+2*nrows) + 2*(ncols-1));
   *mul += nrows*(2 + 4*(ncols-1)*ncols + 2*(ncols-1));
}

void flopcount_dbl_beta_times_V ( int nrows, long long int *mul )
{
   const int nbthreads = nrows;
   // equals rowdim as in the call to the dbl_beta_times kernel

   *mul += nbthreads;
}

void flopcount_cmplx_beta_times_V ( int nrows, long long int *mul )
{
   const int nbthreads = nrows;
   // equals rowdim as in the call to the dbl_beta_times kernel

   *mul += 2*nbthreads;
}

void flopcount_dbl_initialize_WYT ( int dim, long long int *mul )
{
   // the number of threads equals dim*dim
   *mul += dim*dim;
}

void flopcount_cmplx_initialize_WYH
 ( int dim, long long int *add, long long int *mul )
{
   // the number of threads equals dim*dim
   *add += 2*dim*dim;
   *mul += 4*dim*dim;
}

void flopcount_dbl_update_WYT
 ( int dim, long long int *add, long long int *mul )
{
   // the number of threads equals dim*dim
   *add += dim*dim;
   *mul += dim*dim;
}

void flopcount_cmplx_update_WYH
 ( int dim, long long int *add, long long int *mul )
{
   // the number of threads equals dim*dim
   *add += 4*dim*dim;
   *mul += 4*dim*dim;
}

void flopcount_dbl_beta_next_W
 ( int nrows, long long int *add, long long int *mul )
{
   // the number of threads equals nrows
   *add += nrows;
   *mul += nrows + 1;
}

void flopcount_cmplx_beta_next_W
 ( int nrows, long long int *add, long long int *mul )
{
   // the number of threads equals nrows
   *add += 4*nrows;
   *mul += 4*nrows + 2;
}

void flopcount_dbl_small_WYT
 ( int nrows, int szt, long long int *add, long long int *mul )
{
   const int nbthreads = nrows*nrows;

   *add += nbthreads;
   *mul += nbthreads;
}

void flopcount_cmplx_small_WYH
 ( int nrows, int szt, long long int *add, long long int *mul )
{
   const int nbthreads = nrows*nrows;

   *add += 4*nbthreads;
   *mul += 4*nbthreads;
}

void flopcount_dbl_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   long long int *add, long long int *mul )
{
   const int nbthreads = dim*rowdim;

   *add += rowdim*nbthreads;
   *mul += rowdim*nbthreads;
}

void flopcount_cmplx_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   long long int *add, long long int *mul )
{
   const int nbthreads = dim*rowdim;

   *add += 4*rowdim*nbthreads;
   *mul += 4*rowdim*nbthreads;
}

void flopcount_dbl_small_YWTC
 ( int rowdim, int coldim, long long int *add, long long int *mul )
{
   const int nbthreads = rowdim*coldim;

   *add += nbthreads;
   *mul += nbthreads;
}

void flopcount_cmplx_small_YWHC
 ( int rowdim, int coldim, long long int *add, long long int *mul )
{
   const int nbthreads = rowdim*coldim;

   *add += 4*nbthreads;
   *mul += 4*nbthreads;
}

void flopcount_dbl_small_Qupdate
 ( int dim, int rowdim, long long int *add )
{
   const int nbthreads = dim*rowdim;

   *add += nbthreads;
}

void flopcount_cmplx_small_Qupdate
 ( int dim, int rowdim, long long int *add )
{
   const int nbthreads = dim*rowdim;

   *add += 2*nbthreads;
}

void flopcount_dbl_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   long long int *add )
{
   const int rowdim = nrows - rowoff;
   const int nbthreads = rowdim*coldim;

   *add += nbthreads;
}

void flopcount_cmplx_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   long long int *add )
{
   const int rowdim = nrows - rowoff;
   const int nbthreads = rowdim*coldim;

   *add += 2*nbthreads;
}
