// The file dbl8_factors_testers.h specifies test functions
// on matrix factorizations in octo double precision.

#ifndef __dbl8_factors_testers_h__
#define __dbl8_factors_testers_h__

void test_factors_real8_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization on real data. */

void test_factors_cmplx8_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization
 *   on complex data. */

int test_real8_qr_factors
 ( int nrows, int ncols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
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
 *   Ahihihi  has the highest doubles of A;
 *   Alohihi  has the second highest doubles of A; 
 *   Ahilohi  has the third highest doubles of A; 
 *   Alolohi  has the fourth highest doubles of A; 
 *   Ahihilo  has the fourth lowest doubles of A;
 *   Alohilo  has the third lowest doubles of A;
 *   Ahilolo  has the second lowest doubles of A;
 *   Alololo  has the lowest doubles of A;
 *   Qhihihi  has the highest doubles of Q;
 *   Qlohihi  has the second highest doubles of Q;
 *   Qhilohi  has the third highest doubles of Q;
 *   Qlolohi  has the fourth highest doubles of Q;
 *   Qhihilo  has the fourth lowest doubles of Q;
 *   Qlohilo  has the third lowest doubles of Q;
 *   Qhilolo  has the second lowest doubles of Q;
 *   Qlololo  has the lowest doubles of Q;
 *   Rhihihi  has the highest doubles of R;
 *   Rlohihi  has the second highest doubles of R;
 *   Rhilohi  has the third highest doubles of R;
 *   Rlolohi  has the fourth highest doubles of R;
 *   Rhihilo  has the fourth lowest doubles of R;
 *   Rlohilo  has the third lowest doubles of R;
 *   Rhilolo  has the second lowest doubles of R;
 *   Rlololo  has the lowest doubles of R;
 *   tol      is the tolerance on the errors;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_real8_qr_factors_probe
 ( int nrows, int ncols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
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
 *   Ahihihi  has the highest doubles of A;
 *   Alohihi  has the second highest doubles of A; 
 *   Ahilohi  has the third highest doubles of A; 
 *   Alolohi  has the fourth highest doubles of A; 
 *   Ahihilo  has the fourth lowest doubles of A;
 *   Alohilo  has the third lowest doubles of A;
 *   Ahilolo  has the second lowest doubles of A;
 *   Alololo  has the lowest doubles of A;
 *   Qhihihi  has the highest doubles of Q;
 *   Qlohihi  has the second highest doubles of Q;
 *   Qhilohi  has the third highest doubles of Q;
 *   Qlolohi  has the fourth highest doubles of Q;
 *   Qhihilo  has the fourth lowest doubles of Q;
 *   Qlohilo  has the third lowest doubles of Q;
 *   Qhilolo  has the second lowest doubles of Q;
 *   Qlololo  has the lowest doubles of Q;
 *   Rhihihi  has the highest doubles of R;
 *   Rlohihi  has the second highest doubles of R;
 *   Rhilohi  has the third highest doubles of R;
 *   Rlolohi  has the fourth highest doubles of R;
 *   Rhihilo  has the fourth lowest doubles of R;
 *   Rlohilo  has the third lowest doubles of R;
 *   Rhilolo  has the second lowest doubles of R;
 *   Rlololo  has the lowest doubles of R;
 *   nbprobes equals the number of probes;
 *   verbose  is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_cmplx8_qr_factors
 ( int nrows, int ncols,
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
   double tol, int verbose );
/*
 * DESCRIPTION :
 *   Computes the errors of |Q^T*A - R| and |Q^T*Q - I|, for complex data.
 *
 * NOTE : for large values of nrows and ncols, this test may become
 *   too time consuming, for large dimensions, use the _probe test.
 *
 * ON ENTRY :
 *   nrows     number of rows of A, of R, of Q,
 *             and the number of columns of Q;
 *   ncols     number of columns in A and of R;
 *   Arehihihi has the highest doubles of the real parts of A;
 *   Arelohihi has the second highest doubles of the real parts of A;
 *   Arehilohi has the third highest doubles of the real parts of A;
 *   Arelolohi has the fourth highest doubles of the real parts of A;
 *   Arehihilo has the fourth lowest doubles of the real parts of A;
 *   Arelohilo has the third lowest doubles of the real parts of A;
 *   Arehilolo has the second lowest doubles of the real parts of A;
 *   Arelololo has the lowest doubles of the real parts of A;
 *   Aimhihihi has the highest doubles of the imaginary parts of A;
 *   Aimlohihi has the second highest doubles of the imaginary parts of A;
 *   Aimhilohi has the third highest doubles of the imaginary parts of A;
 *   Aimlolohi has the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo has the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo has the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo has the second lowest doubles of the imaginary parts of A;
 *   Aimlololo has the lowest doubles of the imaginary parts of A;
 *   Qrehihihi has the highest doubles of the real parts of Q;
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
 *   Rrehihihi has the highest doubles of the real parts of R;
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
 *   tol       is the tolerance on the errors;
 *   verbose   is the verbose level.
 *
 * RETURN :
 *   0        if all errors are less than the tolerance tol;
 *   1        if one error is larger than the tolearnce tol. */

int test_cmplx8_qr_factors_probe
 ( int nrows, int ncols,
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
 *   nrows     number of rows in A, in R, and in Q;
 *   ncols     number of columns in A and in R;
 *   Arehihihi has the highest doubles of the real parts of A;
 *   Arelohihi has the second highest doubles of the real parts of A;
 *   Arehilohi has the third highest doubles of the real parts of A;
 *   Arelolohi has the fourth highest doubles of the real parts of A;
 *   Arehihilo has the fourth lowest doubles of the real parts of A;
 *   Arelohilo has the third lowest doubles of the real parts of A;
 *   Arehilolo has the second lowest doubles of the real parts of A;
 *   Arelololo has the lowest doubles of the real parts of A;
 *   Aimhihihi has the highest doubles of the imaginary parts of A;
 *   Aimlohihi has the second highest doubles of the imaginary parts of A;
 *   Aimhilohi has the third highest doubles of the imaginary parts of A;
 *   Aimlolohi has the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo has the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo has the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo has the second lowest doubles of the imaginary parts of A;
 *   Aimlololo has the lowest doubles of the imaginary parts of A;
 *   Qrehihihi has the highest doubles of the real parts of Q;
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
 *   Rrehihihi has the highest doubles of the real parts of R;
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
 *   tol       is the tolerance on the errors;
 *   nbprobes  equals the number of probes;
 *   verbose   is the verbose level.
 *
 * RETURN :
 *   0         if all errors are less than the tolerance tol;
 *   1         if one error is larger than the tolearnce tol. */

void test_factors_real8_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on real data. */

void test_factors_cmplx8_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on complex data. */

#endif
