// The file dbl_factors_testers.h specifies test functions
// on matrix factorizations in double precision.

#ifndef __dbl_factors_testers_h__
#define __dbl_factors_testers_h__

void test_factors_real_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization on real data. */

void test_factors_cmplx_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization
 *   on complex data. */

int test_real_qr_factors
 ( int nrows, int ncols, double **A, double **Q, double **R,
   double tol, int verbose );
/*
 * DESCRIPTION :
 *   Computes the errors of |Q^T*A - R| and |Q^T*Q - I|, for real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in A, in R, and in Q;
 *   ncols    number of columns in A and in R;
 *   A        an nrows-by-ncols matrix;
 *   Q        an nrows-by-nrows matrix;
 *   R        an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_real_qr_factors_probe
 ( int nrows, int ncols, double **A, double **Q, double **R,
   double tol, int nbprobes, int verbose );
/*
 * DESCRIPTION :
 *   Computes the errors of |Q^T*A - R| and |Q^T*Q - I|,
 *   for some random indices, on real data.
 *
 * NOTE : this is more economical for large dimensions,
 *   than the complete comparisons.
 *
 * ON ENTRY :
 *   nrows    number of rows in A, in R, and in Q;
 *   ncols    number of columns in A and in R;
 *   A        an nrows-by-ncols matrix;
 *   Q        an nrows-by-nrows matrix;
 *   R        an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   nbprobes equals the number of probes;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_cmplx_qr_factors
 ( int nrows, int ncols, double **Are, double **Aim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double tol, int verbose );
/*
 * DESCRIPTION :
 *   Computes the errors of |Q^T*A - R| and |Q^T*Q - I|, for complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in A, in R, and in Q;
 *   ncols    number of columns in A and in R;
 *   Are      real parts of an nrows-by-ncols matrix;
 *   Aim      imaginary parts of an nrows-by-ncols matrix;
 *   Qre      real parts of an nrows-by-nrows matrix;
 *   Qim      imaginary parts of an nrows-by-nrows matrix;
 *   Rre      real parts of an nrows-by-ncols matrix;
 *   Rim      imaginary parts of an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_cmplx_qr_factors_probe
 ( int nrows, int ncols, double **Are, double **Aim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double tol, int nbprobes, int verbose );
/*
 * DESCRIPTION :
 *   Computes the errors of |Q^T*A - R| and |Q^T*Q - I|,
 *   for some random indices, on complex data,
 *
 * NOTE : this is more economical for large dimensions,
 *   than the complete comparisons.
 *
 * ON ENTRY :
 *   nrows    number of rows in A, in R, and in Q;
 *   ncols    number of columns in A and in R;
 *   Are      real parts of an nrows-by-ncols matrix;
 *   Aim      imaginary parts of an nrows-by-ncols matrix;
 *   Qre      real parts of an nrows-by-nrows matrix;
 *   Qim      imaginary parts of an nrows-by-nrows matrix;
 *   Rre      real parts of an nrows-by-ncols matrix;
 *   Rim      imaginary parts of an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   nbprobes equals the number of probes;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

void test_factors_real_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on real data. */

void test_factors_cmplx_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on complex data. */

#endif
