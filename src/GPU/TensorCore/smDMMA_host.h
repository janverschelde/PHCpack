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

void random_double_double_matrices
 ( double *A, double *B, double *C,
   int numArows, int numAcols, int numBrows, int numBcols,
   int numCrows, int numCcols, int vrblvl=0 );
/*
 * Defines single indexed matrices which represent rewritten
 * double double matrices.
 *
 * ON ENTRY :
 *   A        memory allocated for numArows rows and numAcols columns;
 *   B        memory allocated for numBrows rows and numBcols columns;
 *   C        memory allocated for numCrows rows and numCcols columns;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   A        rewritten (numArows/12)-by-(numAcols/12) matrix,
 *            padded with some zero rows and columns;
 *   B        rewritten (numBrows/12)-by-numBcols matrix,
 *            padded with some zero rows;
 *   C        initialized all numCrows x numCcols to zero. */

void matMultiplyOnHost
 ( double *A, double *B, double *C, float alpha, float beta,
   int numARows, int numAColumns, int numBRows, int numBColumns,
   int numCRows, int numCColumns );
/*
 * Multiplies the matrices A and B into C. */

#endif
