/* Prototypes of the kernels for a double double matrix multiplication. */

#ifndef __DDMM_KERNELS_H__
#define __DDMM_KERNELS_H__

__global__ void ddmm
 ( int nrows, int ncols, int dim,
   double *Ahi, double *Alo, double *Bhi, double *Blo,
   double *Chi, double *Clo );
/*
 * Given in (Ahi, Alo) is an nrows-by-dim row major matrix A, and
 * given in (Bhi, Blo) is an ncols-by-dim column major matrix B,
 * returns in (Chi, Clo) the product of A with B.
 * All matrices are stored as single indexed arrays. */

void GPU_dd_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **BThi, double **BTlo,
   double **Chi, double **Clo, float* milliseconds );
/*
 * Given in (Ahi, Alo) is an nrows-by-dim row major matrix A, and
 * given in (BThi, BTlo) is an ncols-by-dim column major matrix BT,
 * launches kernels to make the product of A with B. 
 * All matrices on entry and on return are double indexed.
 * Returns in milliseconds the kernel time. */

#endif
