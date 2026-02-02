/* Code to initialize and multiply the test matrices on the host. */

#ifndef __SMDMMA_HOST_H__
#define __SMDMMA_HOST_H__

void init_host_matrices ( double *a, double *b, double *c, int nbrange );
/*
 * Using the definitions in smDMMA_dims.h,
 * defines the matrices a, b, and c, as follows:
 *   a is upper triangular with 1 on diagonal and above the diagonal,
 *     row major with M_GLOBAL rows and K_GLOBAL columns,
 *   b is upper triangular with increasing sequences of nbrange numbers,
 *     each time doubling the previous number,
 *     column major with K_GLOBAL rows and N_GLOBAL columns,
 *   c is the zero matrix, M_GLOBAL rows, N_GLOBAL columns. */

void matMultiplyOnHost
 ( double *A, double *B, double *C, float alpha, float beta,
   int numARows, int numAColumns, int numBRows, int numBColumns,
   int numCRows, int numCColumns );
/*
 * Multiplies the matrices A and B into C. */

#endif
