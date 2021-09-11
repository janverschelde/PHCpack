/* The file dbl8_factorizations.h specifies functions to factor matrices
 * in octo double precision. */

#ifndef __dbl8_factorizations_h__
#define __dbl8_factorizations_h__

void CPU_dbl8_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Bhihihi, double **Blohihi, double **Bhilohi, double **Blolohi,
   double **Bhihilo, double **Blohilo, double **Bhilolo, double **Blololo,
   double **Chihihi, double **Clohihi, double **Chilohi, double **Clolohi,
   double **Chihilo, double **Clohilo, double **Chilolo, double **Clololo );
/*
 * DESCRIPTION :
 *   Computes the product C of the matrix A with B on real data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices A and C;
 *   dim      the number of columns in A and rows in B;
 *   cols     the number of columns in the matrices B and C;
 *   Ahihihi  are the highest doubles of the matrix A;
 *   Alohihi  are the second highest doubles of A;
 *   Ahilohi  are the third highest doubles of A;
 *   Alolohi  are the fourth highest doubles of A;
 *   Ahihilo  are the fourth lowest doubles of A;
 *   Alohilo  are the third lowest doubles of A;
 *   Ahilolo  are the second lowest doubles of A;
 *   Alololo  are the lowest doubles of A;
 *   Bhihihi  are the highest doubles of the matrix B;
 *   Blohihi  are the second highest doubles of B;
 *   Bhilohi  are the third highest doubles of B;
 *   Blolohi  are the fourth highest doubles of B;
 *   Bhihilo  are the fourth lowest doubles of B;
 *   Blohilo  are the third lowest doubles of B;
 *   Bhilolo  are the second lowest doubles of B;
 *   Blololo  are the lowest doubles of B;
 *   Chihihi  has space allocated for a rows-by-cols matrix;
 *   Clohihi  has space allocated for a rows-by-cols matrix;
 *   Chilohi  has space allocated for a rows-by-cols matrix.
 *   Clolohi  has space allocated for a rows-by-cols matrix;
 *   Chihilo  has space allocated for a rows-by-cols matrix;
 *   Clohilo  has space allocated for a rows-by-cols matrix;
 *   Chilolo  has space allocated for a rows-by-cols matrix;
 *   Clololo  has space allocated for a rows-by-cols matrix.
 *
 * ON RETURN :
 *   Chihihi  are the highest doubles of C, the product of A with B;
 *   Clohihi  are the second highest doubles of C;
 *   Chilohi  are the third highest doubles of C;
 *   Clolohi  are the fourth highest doubles of C;
 *   Chihilo  are the fourth lowest doubles of C;
 *   Clohilo  are the third lowest doubles of C;
 *   Chilolo  are the second lowest doubles of C;
 *   Clololo  are the lowest doubles of C. */

void CPU_cmplx8_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Brehihihi, double **Brelohihi,
   double **Brehilohi, double **Brelolohi,
   double **Brehihilo, double **Brelohilo,
   double **Brehilolo, double **Brelololo,
   double **Bimhihihi, double **Bimlohihi,
   double **Bimhilohi, double **Bimlolohi,
   double **Bimhihilo, double **Bimlohilo,
   double **Bimhilolo, double **Bimlololo,
   double **Crehihihi, double **Crelohihi,
   double **Crehilohi, double **Crelolohi,
   double **Crehihilo, double **Crelohilo,
   double **Crehilolo, double **Crelololo,
   double **Cimhihihi, double **Cimlohihi,
   double **Cimhilohi, double **Cimlolohi,
   double **Cimhihilo, double **Cimlohilo,
   double **Cimhilolo, double **Cimlololo );
/*
 * DESCRIPTION :
 *   Computes the product C of the matrix A with B on complex data.
 *
 * ON ENTRY :
 *   rows      the number of rows in the matrices A and C;
 *   dim       the number of columns in A and rows in B;
 *   cols      the number of columns in the matrices B and C;
 *   Arehihihi are the highest doubles of the real parts of A;
 *   Arelohihi are the second highest doubles of the real parts of A;
 *   Arehilohi are the third highest doubles of the real parts of A;
 *   Arelolohi are the fourth highest doubles of the real parts of A;
 *   Arehihilo are the fourth lowest doubles of the real parts of A;
 *   Arelohilo are the third lowest doubles of the real parts of A;
 *   Arehilolo are the second lowest doubles of the real parts of A;
 *   Arelololo are the lowest doubles of the real parts of A;
 *   Aimhihihi are the highest doubles of the imaginary parts of A; 
 *   Aimlohihi are the second highest doubles of the imaginary parts of A; 
 *   Aimhilohi are the third highest doubles of the imaginary parts of A; 
 *   Aimlolohi are the fourth highest doubles of the imaginary parts of A; 
 *   Aimhihilo are the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo are the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo are the second lowest doubles of the imaginary parts of A;
 *   Aimlololo are the lowest doubles of the imaginary parts of A;
 *   Brehihihi are the highest doubles of the real parts of B;
 *   Brelohihi are the second highest doubles of the real parts of B;
 *   Brehilohi are the thirdd highest doubles of the real parts of B;
 *   Brelolohi are the fourth highest doubles of the real parts of B;
 *   Brehihilo are the fourth lowest doubles of the real parts of B;
 *   Brelohilo are the third lowest doubles of the real parts of B;
 *   Brehilolo are the second lowest doubles of the real parts of B;
 *   Brelololo are the lowest doubles of the real parts of B;
 *   Bimhihihi are the highest doubles of the imaginary parts of B;
 *   Bimlohihi are the second highest doubles of the imaginary parts of B;
 *   Bimhilohi are the third highest doubles of the imaginary parts of B;
 *   Bimlolohi are the fourth highest doubles of the imaginary parts of B;
 *   Bimhihilo are the fourth lowest doubles of the imaginary parts of B;
 *   Bimlohilo are the third lowest doubles of the imaginary parts of B;
 *   Bimhilolo are the second lowest doubles of the imaginary parts of B;
 *   Bimlololo are the lowest doubles of the imaginary parts of B;
 *   Crehihihi has space allocated for a rows-by-cols matrix;
 *   Crelohihi has space allocated for a rows-by-cols matrix;
 *   Crehilohi has space allocated for a rows-by-cols matrix;
 *   Crelolohi has space allocated for a rows-by-cols matrix;
 *   Crehihilo has space allocated for a rows-by-cols matrix;
 *   Crelohilo has space allocated for a rows-by-cols matrix;
 *   Crehilolo has space allocated for a rows-by-cols matrix;
 *   Crelololo has space allocated for a rows-by-cols matrix;
 *   Cimhihihi has space allocated for a rows-by-cols matrix;
 *   Cimlohihi has space allocated for a rows-by-cols matrix;
 *   Cimhilohi has space allocated for a rows-by-cols matrix;
 *   Cimlolohi has space allocated for a rows-by-cols matrix.
 *   Cimhihilo has space allocated for a rows-by-cols matrix;
 *   Cimlohilo has space allocated for a rows-by-cols matrix;
 *   Cimhilolo has space allocated for a rows-by-cols matrix;
 *   Cimlololo has space allocated for a rows-by-cols matrix.
 *
 * ON RETURN :
 *   Crehihihi are the highest doubles of the real parts of C,
               the product of A with B;
 *   Crelohihi are the second highest doubles of the real parts of C;
 *   Crehilohi are the third highest doubles of the real parts of C;
 *   Crelolohi are the fourth highest doubles of the real parts of C;
 *   Crehihilo are the fourth lowest doubles of the real parts of C;
 *   Crelohilo are the thirdd lowest doubles of the real parts of C;
 *   Crehilolo are the second lowest doubles of the real parts of C;
 *   Crelololo are the lowest doubles of the real parts of C;
 *   Cimhihihi are the highest doubles of the imaginary parts of C;
 *   Cimlohihi are the second highest doubles of the imaginary parts of C;
 *   Cimhilohi are the third highest doubles of the imaginary parts of C;
 *   Cimlolohi are the fourth highest doubles of the imaginary parts of C;
 *   Cimhihilo are the fourth lowest doubles of the imaginary parts of C;
 *   Cimlohilo are the third lowest doubles of the imaginary parts of C;
 *   Cimhilolo are the second lowest doubles of the imaginary parts of C;
 *   Cimlololo are the lowest doubles of the imaginary parts of C. */

void CPU_dbl8_factors_forward
 ( int dim,
   double **Lhihihi, double **Llohihi, double **Lhilohi, double **Llolohi,
   double **Lhihilo, double **Llohilo, double **Lhilolo, double **Llololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo );
/*
 * DESCRIPTION :
 *   Solves the lower triangular system L*x = b with forward substitution,
 *   on real data.
 *
 * REQUIRED : the matrix L has ones on the diagonal.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix L;
 *   Lhihihi  has the highest doubles of L;
 *   Llohihi  has the second highest doubles of L;
 *   Llohihi  has the third highest doubles of L;
 *   Llohihi  has the fourth highest doubles of L;
 *   Lhilolo  has the fourth lowest doubles of L;
 *   Lhilolo  has the third lowest doubles of L;
 *   Lhilolo  has the second lowest doubles of L;
 *   Llololo  has the lowest doubles of L;
 *   bhihihi  has the highest doubles of b;
 *   blohihi  has the second highest doubles of b;
 *   blohihi  has the third highest doubles of b;
 *   blohihi  has the fourth highest doubles of b;
 *   bhilolo  has the fourth lowest doubles of b;
 *   bhilolo  has the third lowest doubles of b;
 *   bhilolo  has the second lowest doubles of b;
 *   blololo  has the lowest doubles of b;
 *   xhihihi  has space for dim doubles;
 *   xlohihi  has space for dim doubles;
 *   xhilohi  has space for dim doubles;
 *   xlolohi  has space for dim doubles;
 *   xhihilo  has space for dim doubles;
 *   xlohilo  has space for dim doubles;
 *   xhilolo  has space for dim doubles;
 *   xlololo  has space for dim doubles.
 *
 * ON RETURN :
 *   xhihihi  has the highest doubles of the solution x;
 *   xlohihi  has the second highest doubles of x;
 *   xhilohi  has the third highest doubles of x;
 *   xlolohi  has the fourth highest doubles of x;
 *   xhihilo  has the fourth lowest doubles of x;
 *   xlohilo  has the third lowest doubles of x;
 *   xhilolo  has the second lowest doubles of x;
 *   xlololo  has the lowest doubles of x. */

void CPU_cmplx8_factors_forward
 ( int dim,
   double **Lrehihihi, double **Lrelohihi,
   double **Lrehilohi, double **Lrelolohi,
   double **Lrehihilo, double **Lrelohilo,
   double **Lrehilolo, double **Lrelololo,
   double **Limhihihi, double **Limlohihi,
   double **Limhilohi, double **Limlolohi,
   double **Limhihilo, double **Limlohilo,
   double **Limhilolo, double **Limlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo );
/*
 * DESCRIPTION :
 *   Solves the lower triangular system L*x = b with forward substitution,
 *   on real data.
 *
 * REQUIRED : the matrix L has ones on the diagonal.
 *
 * ON ENTRY :
 *   dim       number of rows and columns in the matrix L;
 *   Lrehihihi has the highest doubles of the real parts of L;
 *   Lrelohihi has the second highest doubles of the real parts of L;
 *   Lrehilohi has the third highest doubles of the real parts of L;
 *   Lrelolohi has the fourth highest doubles of the real parts of L;
 *   Lrehihilo has the fourth lowest doubles of the real parts of L;
 *   Lrelohilo has the third lowest doubles of the real parts of L;
 *   Lrehilolo has the second lowest doubles of the real parts of L;
 *   Lrelololo has the lowest doubles of the real parts of L;
 *   Limhihihi has the highest doubles of the imaginary parts of L;
 *   Limlohihi has the second highest doubles of the imaginary parts of L;
 *   Limhilohi has the third highest doubles of the imaginary parts of L;
 *   Limlolohi has the fourth highest doubles of the imaginary parts of L;
 *   Limhihilo has the fourth lowest doubles of the imaginary parts of L;
 *   Limlohilo has the third lowest doubles of the imaginary parts of L;
 *   Limhilolo has the second lowest doubles of the imaginary parts of L;
 *   Limlololo has the lowest doubles of the imaginary parts of L;
 *   brehihihi has the highest doubles of the real parts of b;
 *   brelohihi has the second highest doubles of the real parts of b;
 *   brehilohi has the third highest doubles of the real parts of b;
 *   brelolohi has the fourth highest doubles of the real parts of b;
 *   brehihilo has the fourth lowest doubles of the real parts of b;
 *   brelohilo has the third lowest doubles of the real parts of b;
 *   brehilolo has the second lowest doubles of the real parts of b;
 *   brelololo has the lowest doubles of the real parts of b;
 *   bimhihihi has the highest doubles of the imaginary parts of b;
 *   bimlohihi has the second highest doubles of the imaginary parts of b;
 *   bimhilohi has the third highest doubles of the imaginary parts of b;
 *   bimlolohi has the fourth highest doubles of the imaginary parts of b;
 *   bimhihilo has the fourth lowest doubles of the imaginary parts of b;
 *   bimlohilo has the third lowest doubles of the imaginary parts of b;
 *   bimhilolo has the second lowest doubles of the imaginary parts of b;
 *   bimlololo has the lowest doubles of the imaginary parts of b;
 *   xrehihihi has space for dim doubles;
 *   xrelohihi has space for dim doubles;
 *   xrehilohi has space for dim doubles;
 *   xrelolohi has space for dim doubles;
 *   ximhihilo has space for dim doubles;
 *   ximlohilo has space for dim doubles;
 *   ximhilolo has space for dim doubles;
 *   ximlololo has space for dim doubles;
 *   xrehihihi has space for dim doubles;
 *   xrelohihi has space for dim doubles;
 *   xrehilohi has space for dim doubles;
 *   xrelolohi has space for dim doubles;
 *   ximhihilo has space for dim doubles;
 *   ximlohilo has space for dim doubles;
 *   ximhilolo has space for dim doubles;
 *   ximlololo has space for dim doubles.
 *
 * ON RETURN :
 *   xrehihihi has the highest doubles of the real parts of the solution x;
 *   xrelohihi has the second highest doubles of the real parts of x;
 *   xrehilohi has the third highest doubles of the real parts of x;
 *   xrelolohi has the fourth highest doubles of the real parts of x;
 *   xrehihilo has the fourth lowest doubles of the real parts of x;
 *   xrelohilo has the third lowest doubles of the real parts of x;
 *   xrehilolo has the second lowest doubles of the real parts of x;
 *   xrelololo has the lowest doubles of the real parts of x;
 *   ximhihihi has the highest doubles of the imaginary parts of x;
 *   ximlohihi has the second highest doubles of the imaginary parts of x;
 *   ximhilohi has the third highest doubles of the imaginary parts of x;
 *   ximlolohi has the fourth highest doubles of the imaginary parts of x;
 *   ximhihilo has the fourth lowest doubles of the imaginary parts of x;
 *   ximlohilo has the third lowest doubles of the imaginary parts of x;
 *   ximhilolo has the second lowest doubles of the imaginary parts of x;
 *   ximlololo has the lowest doubles of the imaginary parts of x. */

void CPU_dbl8_factors_backward
 ( int dim,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = b with back substitution,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix U;
 *   Uhihihi  are the highest doubles of U
 *   Ulohihi  are the second highest doubles of U;
 *   Uhilohi  are the third highest doubles of U;
 *   Ulolohi  are the fourth highest doubles of U;
 *   Uhihilo  are the fourth lowest doubles of U;
 *   Ulohilo  are the third lowest doubles of U;
 *   Uhilolo  are the second lowest doubles of U;
 *   Ulololo  are the lowest doubles of U;
 *   bhihihi  are the highest doubles of b;
 *   blohihi  are the second highest doubles of b;
 *   bhilohi  are the third highest doubles of b;
 *   blolohi  are the fourth highest doubles of b;
 *   bhihilo  are the fourth lowest doubles of b;
 *   blohilo  are the third lowest doubles of b;
 *   bhilolo  are the second lowest doubles of b;
 *   blololo  are the lowest doubles of b;
 *   xhihihi  space for dim doubles;
 *   xlohihi  space for dim doubles;
 *   xhilohi  space for dim doubles;
 *   xlolohi  space for dim doubles;
 *   xhihilo  space for dim doubles;
 *   xlohilo  space for dim doubles;
 *   xhilolo  space for dim doubles;
 *   xlololo  space for dim doubles.
 *
 * ON RETURN :
 *   xhihihi  are the highest doubles of the solution x;
 *   xlohihi  are the second highest doubles of x;
 *   xhilohi  are the third highest doubles of x;
 *   xlolohi  are the fourth highest doubles of x;
 *   xhihilo  are the fourth lowest doubles of x;
 *   xlohilo  are the third lowest doubles of x;
 *   xhilolo  are the second lowest doubles of x;
 *   xlololo  are the lowest doubles of x. */

void CPU_cmplx8_factors_backward
 ( int dim,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi,
   double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo,
   double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi,
   double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo,
   double *ximhilolo, double *ximlololo );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = b with back substitution,
 *   on real data.
 *
 * ON ENTRY :
 *   dim       number of rows and columns in the upper triangular matrix U;
 *   Urehihkhi are the highest doubles of the real parts of U;
 *   Urelohihi are the second highest doubles of the real parts of U;
 *   Urehilohi are the third highest doubles of the real parts of U;
 *   Urelolohi are the fourth highest doubles of the real parts of U;
 *   Urehihilo are the fourth lowest doubles of the real parts of U;
 *   Urelohilo are the third lowest doubles of the real parts of U;
 *   Urehilolo are the second lowest doubles of the real parts of U;
 *   Urelololo are the lowest doubles of the real parts of U;
 *   Uimhihihi are the highest doubles of the imaginary parts of U;
 *   Uimlohihi are the second highest doubles of the imaginary parts of U;
 *   Uimhilohi are the third highest doubles of the imaginary parts of U;
 *   Uimlolohi are the fourth highest doubles of the imaginary parts of U;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts of U;
 *   Uimlohilo are the third lowest doubles of the imaginary parts of U;
 *   Uimhilolo are the second lowest doubles of the imaginary parts of U;
 *   Uimlololo are the lowest doubles of the imaginary parts of U;
 *   brehihihi has the highest doubles of the real parts of b;
 *   brelohihi has the second highest doubles of the real parts of b;
 *   brehilohi has the third highest doubles of the real parts of b;
 *   brelolohi has the fourth highest doubles of the real parts of b;
 *   brehihilo has the fourth lowest doubles of the real parts of b;
 *   brelohilo has the third lowest doubles of the real parts of b;
 *   brehilolo has the second lowest doubles of the real parts of b;
 *   brelololo has the lowest doubles of the real parts of b;
 *   bimhihihi has the highest doubles of the imaginary parts of b;
 *   bimlohihi has the second highest doubles of the imaginary parts of b;
 *   bimhilohi has the third highest doubles of the imaginary parts of b;
 *   bimlolohi has the fourth highest doubles of the imaginary parts of b;
 *   bimhihilo has the fourth lowest doubles of the imaginary parts of b;
 *   bimlohilo has the third lowest doubles of the imaginary parts of b;
 *   bimhilolo has the second lowest doubles of the imaginary parts of b;
 *   bimlololo has the lowest doubles of the imaginary parts of b;
 *   xrehihihi has space for dim doubles;
 *   xrelohihi has space for dim doubles;
 *   xrehilohi has space for dim doubles;
 *   xrelolohi has space for dim doubles;
 *   ximhihilo has space for dim doubles;
 *   ximlohilo has space for dim doubles;
 *   ximhilolo has space for dim doubles;
 *   ximlololo has space for dim doubles;
 *   xrehihihi has space for dim doubles;
 *   xrelohihi has space for dim doubles;
 *   xrehilohi has space for dim doubles;
 *   xrelolohi has space for dim doubles;
 *   ximhihilo has space for dim doubles;
 *   ximlohilo has space for dim doubles;
 *   ximhilolo has space for dim doubles;
 *   ximlololo has space for dim doubles.
 *
 * ON RETURN :
 *   xrehihihi has the highest doubles of the real parts of the solution x;
 *   xrelohihi has the second highest doubles of the real parts of x;
 *   xrehilohi has the third highest doubles of the real parts of x;
 *   xrelolohi has the fourth highest doubles of the real parts of x;
 *   xrehihilo has the fourth lowest doubles of the real parts of x;
 *   xrelohilo has the third lowest doubles of the real parts of x;
 *   xrehilolo has the second lowest doubles of the real parts of x;
 *   xrelololo has the lowest doubles of the real parts of x;
 *   ximhihihi has the highest doubles of the imaginary parts of x;
 *   ximlohihi has the second highest doubles of the imaginary parts of x;
 *   ximhilohi has the third highest doubles of the imaginary parts of x;
 *   ximlolohi has the fourth highest doubles of the imaginary parts of x;
 *   ximhihilo has the fourth lowest doubles of the imaginary parts of x;
 *   ximlohilo has the third lowest doubles of the imaginary parts of x;
 *   ximhilolo has the second lowest doubles of the imaginary parts of x;
 *   ximlololo has the lowest doubles of the imaginary parts of x. */

void CPU_dbl8_factors_lufac
 ( int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   int *pivots );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Ahihihi  are the highest doubles of A;
 *   Alohihi  are the second highest doubles of A;
 *   Ahilohi  are the third highest doubles of A;
 *   Alolohi  are the fourth highest doubles of A;
 *   Ahihilo  are the fourth lowest doubles of A;
 *   Alohilo  are the third lowest doubles of A;
 *   Ahilolo  are the second lowest doubles of A;
 *   Alololo  are the lowest doubles of A;
 *   pivots   space for dim pivots.
 *
 * ON RETURN :
 *   Ahihihi  the lower triangular part of Ahihihi contains the highest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahihihi has the highest doubles of the row reduced A;
 *   Alohihi  the lower triangular part of Alohihi contains the second highest
 *            doubles of the multipliers and the upper triangular part
 *            of Alohihi has the second highest doubles of the row reduced A;
 *   Ahilohi  the lower triangular part of Ahilohi contains the third highest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahilohi has the third highest doubles of the row reduced A;
 *   Alolohi  the lower triangular part of Alolohi contains the fourth highest
 *            doubles of the multipliers and the upper triangular part
 *            of Alolohi has the fourth highest doubles of the row reduced A;
 *   Ahihilo  the lower triangular part of Ahihilo contains the fourth lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahihilo has the fourth lowest doubles of the row reduced A;
 *   Alohilo  the lower triangular part of Alohilo contains the third lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Alohilo has the third lowest doubles of the row reduced A;
 *   Ahilolo  the lower triangular part of Ahilolo contains the second lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahilolo has the second lowest doubles of the row reduced A;
 *   Alololo  the lower triangular part of Alololo contains the lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Alololo has the lowest doubles of the row reduced A;
 *   pivots   are the pivots used. */

void CPU_cmplx8_factors_lufac
 ( int dim,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   int *pivots );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   on real data.
 *
 * ON ENTRY :
 *   dim       number of rows and columns in the matrix A;
 *   Arehihihi are the highest doubles of the real parts of A;
 *   Arelohihi are the second highest doubles of the real parts of A;
 *   Arehilohi are the third highest doubles of the real parts of A;
 *   Arelolohi are the fourth highest doubles of the real parts of A;
 *   Arehihilo are the fourth lowest doubles of the real parts of A;
 *   Arelohilo are the third lowest doubles of the real parts of A;
 *   Arehilolo are the second lowest doubles of the real parts of A;
 *   Arelololo are the lowest doubles of the real parts of A;
 *   Aimhihihi are the highest doubles of the imaginary parts of A;
 *   Aimlohihi are the second highest doubles of the imaginary parts of A;
 *   Aimhilohi are the third highest doubles of the imaginary parts of A;
 *   Aimlolohi are the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo are the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo are the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo are the second lowest doubles of the imaginary parts of A;
 *   Aimlololo are the lowest doubles of the imaginary parts of A;
 *   pivots   space for dim pivots.
 *
 * ON RETURN :
 *   Arehihihi are the highest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arelohihi are the second highest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arehilohi are the third highest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arelolohi are the fourth highest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arehihilo are the fourth lowest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arelohilo are the third lowest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arehilolo are the second lowest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arelololo are the lowest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Aimhihihi are the highest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimlohihi are the second highest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimhilohi are the third highest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimlolohi are the fourth highest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimhihilo are the fourth lowest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimlohilo are the third lowest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimhilolo are the second lowest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimlololo are the lowest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   pivots    are the pivots used. */

void CPU_dbl8_factors_lusolve
 ( int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   int *pivots,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Ahihihi  are the highest doubles of A;
 *   Alohihi  are the second highest doubles of A;
 *   Ahilohi  are the third highest doubles of A;
 *   Alolohi  are the fourth highest doubles of A;
 *   Ahihilo  are the fourth lowest doubles of A;
 *   Alohilo  are the third lowest doubles of A;
 *   Ahilolo  are the second lowest doubles of A;
 *   Alololo  are the lowest doubles of A;
 *   pivots   space for dim pivots;
 *   blohihi  are the second highest doubles of b;
 *   bhilohi  are the third highest doubles of b;
 *   blolohi  are the fourth highest doubles of b;
 *   bhihilo  are the fourth lowest doubles of b;
 *   blohilo  are the third lowest doubles of b;
 *   bhilolo  are the second lowest doubles of b;
 *   blololo  are the lowest doubles of b;
 *   xhihihi  space for dim doubles;
 *   xlohihi  space for dim doubles;
 *   xhilohi  space for dim doubles;
 *   xlolohi  space for dim doubles;
 *   xhihilo  space for dim doubles;
 *   xlohilo  space for dim doubles;
 *   xhilolo  space for dim doubles;
 *   xlololo  space for dim doubles.
 *
 * ON RETURN :
 *   Ahihihi  the lower triangular part of Ahihihi contains the highest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahihihi has the highest doubles of the row reduced A;
 *   Alohihi  the lower triangular part of Alohihi contains the second highest
 *            doubles of the multipliers and the upper triangular part
 *            of Alohihi has the second highest doubles of the row reduced A;
 *   Ahilohi  the lower triangular part of Ahilohi contains the third highest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahilohi has the third highest doubles of the row reduced A;
 *   Alolohi  the lower triangular part of Alolohi contains the fourth highest
 *            doubles of the multipliers and the upper triangular part
 *            of Alolohi has the fourth highest doubles of the row reduced A;
 *   Ahihilo  the lower triangular part of Ahihilo contains the fourth lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahihilo has the fourth lowest doubles of the row reduced A;
 *   Alohilo  the lower triangular part of Alohilo contains the third lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Alohilo has the third lowest doubles of the row reduced A;
 *   Ahilolo  the lower triangular part of Ahilolo contains the second lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahilolo has the second lowest doubles of the row reduced A;
 *   Alololo  the lower triangular part of Alololo contains the lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Alololo has the lowest doubles of the row reduced A;
 *   bhihihi  used as work space;
 *   blohihi  used as work space;
 *   bhilohi  used as work space;
 *   blolohi  used as work space;
 *   bhihilo  used as work space;
 *   blohilo  used as work space;
 *   bhilolo  used as work space;
 *   blololo  used as work space;
 *   xhihihi  are the highest doubles of the solution x;
 *   xlohihi  are the second highest doubles of x;
 *   xhilohi  are the third highest doubles of x;
 *   xlolohi  are the fourth highest doubles of x;
 *   xhihilo  are the fourth lowest doubles of x;
 *   xlohilo  are the third lowest doubles of x;
 *   xhilolo  are the second lowest doubles of x;
 *   xlololo  are the lowest doubles of x. */

void CPU_cmplx8_factors_lusolve
 ( int dim,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   int *pivots,
   double *brehihihi, double *brelohihi,
   double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo,
   double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi,
   double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo,
   double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi,
   double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo,
   double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi,
   double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo,
   double *ximhilolo, double *ximlololo );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b, on complex data.
 *
 * ON ENTRY :
 *   dim       number of rows and columns in the matrix A;
 *   Arehihihi are the highest doubles of the real parts of A;
 *   Arelohihi are the second highest doubles of the real parts of A;
 *   Arehilohi are the third highest doubles of the real parts of A;
 *   Arelolohi are the fourth highest doubles of the real parts of A;
 *   Arehihilo are the fourth lowest doubles of the real parts of A;
 *   Arelohilo are the third lowest doubles of the real parts of A;
 *   Arehilolo are the second lowest doubles of the real parts of A;
 *   Arelololo are the lowest doubles of the real parts of A;
 *   Aimhihihi are the highest doubles of the imaginary parts of A;
 *   Aimlohihi are the second highest doubles of the imaginary parts of A;
 *   Aimhilohi are the third highest doubles of the imaginary parts of A;
 *   Aimlolohi are the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo are the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo are the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo are the second lowest doubles of the imaginary parts of A;
 *   Aimlololo are the lowest doubles of the imaginary parts of A;
 *   pivots    space for dim pivots;
 *   brehihihi has the highest doubles of the real parts of b;
 *   brelohihi has the second highest doubles of the real parts of b;
 *   brehilohi has the third highest doubles of the real parts of b;
 *   brelolohi has the fourth highest doubles of the real parts of b;
 *   brehihilo has the fourth lowest doubles of the real parts of b;
 *   brelohilo has the third lowest doubles of the real parts of b;
 *   brehilolo has the second lowest doubles of the real parts of b;
 *   brelololo has the lowest doubles of the real parts of b;
 *   bimhihihi has the highest doubles of the imaginary parts of b;
 *   bimlohihi has the second highest doubles of the imaginary parts of b;
 *   bimhilohi has the third highest doubles of the imaginary parts of b;
 *   bimlolohi has the fourth highest doubles of the imaginary parts of b;
 *   bimhihilo has the fourth lowest doubles of the imaginary parts of b;
 *   bimlohilo has the third lowest doubles of the imaginary parts of b;
 *   bimhilolo has the second lowest doubles of the imaginary parts of b;
 *   bimlololo has the lowest doubles of the imaginary parts of b;
 *   xrehihihi has space for dim doubles;
 *   xrelohihi has space for dim doubles;
 *   xrehilohi has space for dim doubles;
 *   xrelolohi has space for dim doubles;
 *   ximhihilo has space for dim doubles;
 *   ximlohilo has space for dim doubles;
 *   ximhilolo has space for dim doubles;
 *   ximlololo has space for dim doubles;
 *   xrehihihi has space for dim doubles;
 *   xrelohihi has space for dim doubles;
 *   xrehilohi has space for dim doubles;
 *   xrelolohi has space for dim doubles;
 *   ximhihilo has space for dim doubles;
 *   ximlohilo has space for dim doubles;
 *   ximhilolo has space for dim doubles;
 *   ximlololo has space for dim doubles.
 *
 * ON RETURN :
 *   Arehihihi are the highest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arelohihi are the second highest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arehilohi are the third highest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arelolohi are the fourth highest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arehihilo are the fourth lowest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arelohilo are the third lowest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arehilolo are the second lowest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Arelololo are the lowest doubles of the real parts
 *             of the multipliers and the row reduced A;
 *   Aimhihihi are the highest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimlohihi are the second highest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimhilohi are the third highest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimlolohi are the fourth highest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimhihilo are the fourth lowest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimlohilo are the third lowest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimhilolo are the second lowest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   Aimlololo are the lowest doubles of the imaginary parts
 *             of the multipliers and the row reduced A;
 *   brehihihi was used as work space;
 *   brelohihi was used as work space;
 *   brehilohi was used as work space;
 *   brelolohi was used as work space;
 *   brehihilo was used as work space;
 *   brelohilo was used as work space;
 *   brehilolo was used as work space;
 *   brelololo was used as work space;
 *   bimhihihi was used as work space;
 *   bimlohihi was used as work space;
 *   bimhilohi was used as work space;
 *   bimlolohi was used as work space;
 *   bimhihilo was used as work space;
 *   bimlohilo was used as work space;
 *   bimhilolo was used as work space;
 *   bimlololo was used as work space;
 *   xrehihihi has the highest doubles of the real parts of the solution x;
 *   xrelohihi has the second highest doubles of the real parts of x;
 *   xrehilohi has the third highest doubles of the real parts of x;
 *   xrelolohi has the fourth highest doubles of the real parts of x;
 *   xrehihilo has the fourth lowest doubles of the real parts of x;
 *   xrelohilo has the third lowest doubles of the real parts of x;
 *   xrehilolo has the second lowest doubles of the real parts of x;
 *   xrelololo has the lowest doubles of the real parts of x;
 *   ximhihihi has the highest doubles of the imaginary parts of x;
 *   ximlohihi has the second highest doubles of the imaginary parts of x;
 *   ximhilohi has the third highest doubles of the imaginary parts of x;
 *   ximlolohi has the fourth highest doubles of the imaginary parts of x;
 *   ximhihilo has the fourth lowest doubles of the imaginary parts of x;
 *   ximlohilo has the third lowest doubles of the imaginary parts of x;
 *   ximhilolo has the second lowest doubles of the imaginary parts of x;
 *   ximlololo has the lowest doubles of the imaginary parts of x. */

/*

void CPU_dbl8_factors_house
 ( int n,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo );
*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   xhihi    the n high doubles of x;
 *   xlohi    the n high doubles of x;
 *   xhilo    the n low doubles of x;
 *   xlolo    the n low doubles of x;
 *   vhihi    space for n doubles;
 *   vlohi    space for n doubles;
 *   vhilo    space for n doubles;
 *   vlolo    space for n doubles.
 *
 * ON RETURN :
 *   vhihi    highest doubles of the Householder vector;
 *   vlohi    second highest doubles of the Householder vector;
 *   vhilo    second lowest doubles of the Householder vector;
 *   vlolo    lowest doubles of the Householder vector;
 *   betahihi equals the highest double of 2/(transpose(v)*v);
 *   betalohi equals the second highest double of 2/(transpose(v)*v);
 *   betahilo equals the second lowest double of 2/(transpose(v)*v);
 *   betalolo equals the lowest double of 2/(transpose(v)*v). *

void CPU_cmplx8_factors_house 
( int n,
  double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
  double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
  double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
  double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
  double *betahihi, double *betalohi, double *betahilo, double *betalolo );
*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   xrehihi  highest doubles of the real parts of the vector x;
 *   xrelohi  second highest doubles of the real parts of the vector x;
 *   xrehilo  second lowest doubles of the real parts of the vector x;
 *   xrelolo  lowest doubles of the real parts of the vector x;
 *   ximhihi  highest doubles of the imaginary parts of the vector x;
 *   ximlohi  second highest doubles of the imaginary parts of the vector x;
 *   ximhilo  second lowest doubles of the imaginary parts of the vector x;
 *   ximlolo  lowest doubles of the imaginary parts of the vector x;
 *   vrehihi  space for n doubles;
 *   vrelohi  space for n doubles;
 *   vrehilo  space for n doubles;
 *   vrelolo  space for n doubles;
 *   vimhihi  space for n doubles;
 *   vimlohi  space for n doubles;
 *   vimhilo  space for n doubles;
 *   vimlolo  space for n doubles.
 *
 * ON RETURN :
 *   vrehihi  highest doubles of the real parts of the Householder vector v;
 *   vrelohi  second highest doubles of the real parts of v;
 *   vrehilo  second lowest doubles of the real parts of v;
 *   vrelolo  lowest doubles of the real parts of v;
 *   vimhihi  highest doubles of the imaginary parts of v;
 *   vimlohi  second highest doubles of the imaginary parts of v;
 *   vimhilo  second lowest doubles of the imaginary parts of v;
 *   vimlolo  lowest doubles of the imaginary parts of v;
 *   betahihi is the highest double of 2/(transpose(v)*v);
 *   betalohi is the second highest double of 2/(transpose(v)*v);
 *   betahilo is the second lowest double of 2/(transpose(v)*v);
 *   betalolo is the lowest double of 2/(transpose(v)*v). *

void CPU_dbl8_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double betahihi, double betalohi, double betahilo, double betalolo );
*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   Rhihi    highest doubles of an nrows-by-ncols matrix;
 *   Rlohi    second highest doubles of an nrows-by-ncols matrix;
 *   Rhilo    second lowest doubles of an nrows-by-ncols matrix;
 *   Rlolo    lowest doubles of an nrows-by-ncols matrix;
 *   vhihi    highest doubles of the Householder vector;
 *   vlohi    second highest doubles of the Householder vector;
 *   vhilo    second lowest doubles of the Householder vector;
 *   vlolo    lowest doubles of the Householder vector;
 *   betahihi is the betahihi computed by CPU_dbl8_factors_house;
 *   betalohi is the betalohi computed by CPU_dbl8_factors_house;
 *   betahilo is the betahilo computed by CPU_dbl8_factors_house;
 *   betalolo is the betalolo computed by CPU_dbl8_factors_house.
 *
 * ON RETURN :
 *   Rhihi    highest doubles of the update with the Householder matrix;
 *   Rlohi    second highest doubles of the update;
 *   Rhilo    second lowest doubles of the update;
 *   Rlolo    lowest doubles of the update with the Householder matrix. *

void CPU_cmplx8_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double betahihi, double betalohi, double betahilo, double betalolo );
*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   Rrehihi  highest doubles of the real parts of an nrows-by-ncols matrix R;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R;
 *   vrehihi  highest doubles of the real parts of the Householder vector v;
 *   vrelohi  second highest doubles of the real parts of v;
 *   vrehilo  second lowest doubles of the real parts of v;
 *   vrelolo  lowest doubles of the real parts of v;
 *   vimhihi  highest doubles of the imaginary parts of v;
 *   vimlohi  second highest doubles of the imaginary parts of v;
 *   vimhilo  second lowest doubles of the imaginary parts of v;
 *   vimlolo  lowest doubles of the imaginary parts of v;
 *   betahihi is the betahihi computed by CPU_cmplx8_factors_house;
 *   betalohi is the betalohi computed by CPU_cmplx8_factors_house;
 *   betahilo is the betahilo computed by CPU_cmplx8_factors_house;
 *   betalolo is the betalolo computed by CPU_cmplx8_factors_house.
 *
 * ON RETURN :
 *   Rrehihi  highest doubles of the real parts of the update
 *            with the Householder matrix;
 *   Rrelohi  second highest doubles of the real parts of the update;
 *   Rrehilo  second lowest doubles of the real parts of the update;
 *   Rrelolo  lowest doubles of the real parts of the update;
 *   Rimhihi  highest doubles of the imaginary parts of the update;
 *   Rimlohi  second highest doubles of the imaginary parts of the update;
 *   Rimhilo  second lowest doubles of the imaginary parts of the update;
 *   Rimlolo  lowest doubles of the imaginary parts of the update. *

void CPU_dbl8_factors_rightQupdate
 ( int n, int k,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double betahihi, double betalohi, double betahilo, double betalolo );
*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Qhihi    highest doubles of an n-by-n matrix;
 *   Qlohi    highest doubles of an n-by-n matrix;
 *   Qhilo    lowest doubles of an n-by-n matrix;
 *   Qlolo    lowest doubles of an n-by-n matrix;
 *   vhihi    highest doubles of the Householder vector;
 *   vlohi    second highest doubles of the Householder vector;
 *   vhilo    second lowest doubles of the Householder vector;
 *   vlolo    lowest doubles of the Householder vector;
 *   betahihi is the betahihi computed by CPU_dbl8_factors_house;
 *   betalohi is the betalohi computed by CPU_dbl8_factors_house;
 *   betahilo is the betahilo computed by CPU_dbl8_factors_house;
 *   betalolo is the betalolo computed by CPU_dbl8_factors_house.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of the update with the Householder matrix;
 *   Qlohi    second highest doubles of the update;
 *   Qhilo    second lowest doubles of the update;
 *   Qlolo    lowest doubles of the update with the Householder matrix. *

void CPU_cmplx8_factors_rightQupdate
 ( int n, int k,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double betahihi, double betalohi, double betahilo, double betalolo );
*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Qrehihi  highest doubles of the real parts of an n-by-n matrix Q;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   vrehihi  highest doubles of the real parts of the Householder vector v;
 *   vrelohi  second highest doubles of the real parts of v;
 *   vrehilo  second lowest doubles of the real parts of v;
 *   vrelolo  lowest doubles of the real parts of v;
 *   vimhihi  highest doubles of the imaginary parts of v;
 *   vimlohi  second highest doubles of the imaginary parts of v;
 *   vimhilo  second lowest doubles of the imaginary parts of v;
 *   vimlolo  lowest doubles of the imaginary parts of v;
 *   betahihi is the betahihi computed by CPU_cmplx8_factors_house;
 *   betalohi is the betalohi computed by CPU_cmplx8_factors_house;
 *   betahilo is the betahilo computed by CPU_cmplx8_factors_house;
 *   betalolo is the betalolo computed by CPU_cmplx8_factors_house.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of the update
 *            with the Householder matrix;
 *   Qrelohi  second highest doubles of the real parts of the update;
 *   Qrehilo  second lowest doubles of the real parts of the update;
 *   Qrelolo  lowest doubles of the real parts of the update;
 *   Qimhihi  highest doubles of the imaginary parts of the update;
 *   Qimlohi  second highest doubles of the imaginary parts of the update;
 *   Qimhilo  second lowest doubles of the imaginary parts of the update;
 *   Qimlolo  lowest doubles of the imaginary parts of the update. *

void CPU_dbl8_factors_houseqr
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo );
*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   Ahihi    highest doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Alohi    second highest doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Ahilo    second lowest doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Alolo    lowest doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Qhihi    space for an nrows-by-nrows matrix;
 *   Qlohi    space for an nrows-by-nrows matrix;
 *   Qhilo    space for an nrows-by-nrows matrix;
 *   Qlolo    space for an nrows-by-nrows matrix;
 *   Rhihi    space for an nrows-by-ncols matrix;
 *   Rlohi    space for an nrows-by-ncols matrix;
 *   Rhilo    space for an nrows-by-ncols matrix;
 *   Rlolo    space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Qlohi    second highest doubles of an orthogonal matrix;
 *   Qhilo    second lowest doubles of an orthogonal matrix;
 *   Qlolo    lowest doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Rhihi    highest doubles of the reduced upper triangular form;
 *   Rlohi    second highest doubles of the reduced upper triangular form;
 *   Rhilo    second lowest doubles of the reduced upper triangular form;
 *   Rlolo    lowest doubles of the reduced upper triangular form. *

void CPU_cmplx8_factors_houseqr
 ( int nrows, int ncols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo );
*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   Arehihi  highest doubles of the real parts of an nrows-by-ncols
 *            matrix A, stored as nrows arrays of ncols numbers;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the real parts of :
 *   Arelolo  lowest doubles of the real parts of A;
 *   Aimhi    highest doubles of the imaginary parts of A;
 *   Aimhi    second highest doubles of the imaginary parts of A;
 *   Aimlo    second lowest doubles of the imaginary parts of A;
 *   Aimlo    lowest doubles of the imaginary parts of A;
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
 *   Rimlolo  space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts
 *            of the orthogonal matrix Q, transpose(Q)*A = R;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   Rrehi    highest doubles of the real parts of R,
 *            the reduced upper triangular form of A;
 *   Rrelohi  second highest doubles of the real parts of R,
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R. *

*/

#endif
