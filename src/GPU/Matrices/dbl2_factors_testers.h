// The file dbl2_factors_testers.h specifies test functions
// on matrix factorizations in double double precision.

#ifndef __dbl2_factors_testers_h__
#define __dbl2_factors_testers_h__

void test_factors_real2_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization on real data. */

void test_factors_cmplx2_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization
 *   on complex data. */

int test_real2_qr_factors
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double tol, int verbose );
/*
 * DESCRIPTION :
 *   Computes the errors of |Q^T*A - R| and |Q^T*Q - I|, for real data.
 *
 * NOTE : for large values of nrows and ncols, this test may become
 *   too time consuming, for large dimensions, use the _probe test.
 *
 * ON ENTRY :
 *   nrows    number of rows in A, in R, and in Q;
 *   ncols    number of columns in A and in R;
 *   Ahi      high doubles of an nrows-by-ncols matrix;
 *   Alo      low doubles of an nrows-by-ncols matrix;
 *   Qhi      high doubles of an nrows-by-nrows matrix;
 *   Qlo      low doubles of an nrows-by-nrows matrix;
 *   Rhi      high doubles of an nrows-by-ncols matrix;
 *   Rlo      low doubles of an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_real2_qr_factors_probe
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
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
 *   Ahi      high doubles of an nrows-by-ncols matrix;
 *   Alo      low doubles of an nrows-by-ncols matrix;
 *   Qhi      high doubles of an nrows-by-nrows matrix;
 *   Qlo      low doubles of an nrows-by-nrows matrix;
 *   Rhi      high doubles of an nrows-by-ncols matrix;
 *   Rlo      low doubles of an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   nbprobes equals the number of probes;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_cmplx2_qr_factors
 ( int nrows, int ncols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double tol, int verbose );
/*
 * DESCRIPTION :
 *   Computes the errors of |Q^T*A - R| and |Q^T*Q - I|, for complex data.
 *
 * NOTE : for large values of nrows and ncols, this test may become
 *   too time consuming, for large dimensions, use the _probe test.
 *
 * ON ENTRY :
 *   nrows    number of rows in A, in R, and in Q;
 *   ncols    number of columns in A and in R;
 *   Arehi    high doubles of the real parts of an nrows-by-ncols matrix;
 *   Arelo    low doubles of the real parts of an nrows-by-ncols matrix;
 *   Aimhi    high doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Aimlo    low doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Qrehi    high doubles of the real parts of an nrows-by-nrows matrix;
 *   Qrelo    low doubles of the real parts of an nrows-by-nrows matrix;
 *   Qimhi    high doubles of the imaginary parts of an nrows-by-nrows matrix;
 *   Qimlo    low doubles of the imaginary parts of an nrows-by-nrows matrix;
 *   Rrehi    high doubles of the real parts of an nrows-by-ncols matrix;
 *   Rrelo    low doubles of the real parts of an nrows-by-ncols matrix;
 *   Rimhi    high doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Rimlo    low doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_cmplx2_qr_factors_probe
 ( int nrows, int ncols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
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
 *   Arehi    high doubles of the real parts of an nrows-by-ncols matrix;
 *   Arelo    low doubles of the real parts of an nrows-by-ncols matrix;
 *   Aimhi    high doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Aimlo    low doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Qrehi    high doubles of the real parts of an nrows-by-nrows matrix;
 *   Qrelo    low doubles of the real parts of an nrows-by-nrows matrix;
 *   Qimhi    high doubles of the imaginary parts of an nrows-by-nrows matrix;
 *   Qimlo    low doubles of the imaginary parts of an nrows-by-nrows matrix;
 *   Rrehi    high doubles of the real parts of an nrows-by-ncols matrix;
 *   Rrelo    low doubles of the real parts of an nrows-by-ncols matrix;
 *   Rimhi    high doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Rimlo    low doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   nbprobes equals the number of probes;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

void test_factors_real2_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on real data. */

void test_factors_cmplx2_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on complex data. */

#endif
