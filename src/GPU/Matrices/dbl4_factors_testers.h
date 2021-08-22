// The file dbl4_factors_testers.h specifies test functions
// on matrix factorizations in quad double precision.

#ifndef __dbl4_factors_testers_h__
#define __dbl4_factors_testers_h__

void test_factors_real4_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization on real data. */

void test_factors_cmplx4_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization
 *   on complex data. */

int test_real4_qr_factors
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
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
 *   Ahihi    highest doubles of an nrows-by-ncols matrix;
 *   Alohi    second highest doubles of an nrows-by-ncols matrix;
 *   Ahilo    second lowest doubles of an nrows-by-ncols matrix;
 *   Alolo    lowest doubles of an nrows-by-ncols matrix;
 *   Qhihi    highest doubles of an nrows-by-nrows matrix;
 *   Qlohi    second highest doubles of an nrows-by-nrows matrix;
 *   Qhilo    second lowest doubles of an nrows-by-nrows matrix;
 *   Qlolo    lowest doubles of an nrows-by-nrows matrix;
 *   Rhihi    highest doubles of an nrows-by-ncols matrix;
 *   Rlohi    second highest doubles of an nrows-by-ncols matrix;
 *   Rhilo    second low doubles of an nrows-by-ncols matrix;
 *   Rlolo    low doubles of an nrows-by-ncols matrix;
 *   tol      is the tolerance on the errors;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_real4_qr_factors_probe
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
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
 *   Ahihi    highest doubles of an nrows-by-ncols matrix;
 *   Alohi    second highest doubles of an nrows-by-ncols matrix;
 *   Ahilo    second lowest doubles of an nrows-by-ncols matrix;
 *   Alolo    lowest doubles of an nrows-by-ncols matrix;
 *   Qhihi    highest doubles of an nrows-by-nrows matrix;
 *   Qlohi    second highest doubles of an nrows-by-nrows matrix;
 *   Qhilo    second lowest doubles of an nrows-by-nrows matrix;
 *   Qlolo    lowest doubles of an nrows-by-nrows matrix;
 *   Rhihi    highest doubles of an nrows-by-ncols matrix;
 *   Rlohi    second highest doubles of an nrows-by-ncols matrix;
 *   Rhilo    second low doubles of an nrows-by-ncols matrix;
 *   Rlolo    low doubles of an nrows-by-ncols matrix;
 *   nbprobes equals the number of probes;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_cmplx4_qr_factors
 ( int nrows, int ncols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double tol, int verbose );
/*
 * DESCRIPTION :
 *   Computes the errors of |Q^T*A - R| and |Q^T*Q - I|, for complex data.
 *
 * NOTE : for large values of nrows and ncols, this test may become
 *   too time consuming, for large dimensions, use the _probe test.
 *
 * ON ENTRY :
 *   nrows    number of rows of A, of R, of Q,
 *            and the number of columns of Q;
 *   ncols    number of columns in A and of R;
 *   Arehihi  highest doubles of the real parts of A;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the real parts of A;
 *   Arelolo  lowest doubles of the real parts of A;
 *   Aimhihi  highest doubles of the imaginary parts of A;
 *   Aimlohi  second highest doubles of the imaginary parts of A;
 *   Aimhilo  second lowest doubles of the imaginary parts of A;
 *   Aimlolo  lowest doubles of the imaginary parts of A;
 *   Qrehihi  highest doubles of the real parts of Q;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   Rrehihi  highest doubles of the real parts of R;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R; 
 *   Rrelolo  lowest doubles of the real parts of R; 
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R;
 *   tol      is the tolerance on the errors;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_cmplx4_qr_factors_probe
 ( int nrows, int ncols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
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
 *   Arehihi  highest doubles of the real parts of A;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the real parts of A;
 *   Arelolo  lowest doubles of the real parts of A;
 *   Aimhihi  highest doubles of the imaginary parts of A;
 *   Aimlohi  second highest doubles of the imaginary parts of A;
 *   Aimhilo  second lowest doubles of the imaginary parts of A;
 *   Aimlolo  lowest doubles of the imaginary parts of A;
 *   Qrehihi  highest doubles of the real parts of Q;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   Rrehihi  highest doubles of the real parts of R;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R; 
 *   Rrelolo  lowest doubles of the real parts of R; 
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R;
 *   tol      is the tolerance on the errors;
 *   nbprobes equals the number of probes;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

void test_factors_real4_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on real data. */

void test_factors_cmplx4_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on complex data. */

#endif
