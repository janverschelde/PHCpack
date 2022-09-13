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

void CPU_dbl8_factors_house
 ( int n,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   xhihihi  are the n highest doubles of x;
 *   xlohihi  are the n second highest doubles of x;
 *   xhilohi  are the n third highest doubles of x;
 *   xlolohi  are the n fourth highest doubles of x;
 *   xhihilo  are the n fourth lowest doubles of x;
 *   xlohilo  are the n third lowest doubles of x;
 *   xhilolo  are the n second lowest doubles of x;
 *   xlololo  are the n lowest doubles of x;
 *   vhihihi  has space for n doubles;
 *   vlohihi  has space for n doubles;
 *   vhilohi  has space for n doubles;
 *   vlolohi  has space for n doubles;
 *   vhihilo  has space for n doubles;
 *   vlohilo  has space for n doubles;
 *   vhilolo  has space for n doubles;
 *   vlololo  has space for n doubles.
 *
 * ON RETURN :
 *   vhihihi  highest doubles of the Householder vector;
 *   vlohihi  second highest doubles of the Householder vector;
 *   vhilohi  third highest doubles of the Householder vector;
 *   vlolohi  fourth highest doubles of the Householder vector;
 *   vhihilo  fourth lowest doubles of the Householder vector;
 *   vlohilo  third lowest doubles of the Householder vector;
 *   vhilolo  second lowest doubles of the Householder vector;
 *   vlololo  lowest doubles of the Householder vector;
 *   betahihihi equals the highest double of 2/(transpose(v)*v);
 *   betalohihi equals the second highest double of 2/(transpose(v)*v);
 *   betahilohi equals the third highest double of 2/(transpose(v)*v);
 *   betalolohi equals the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo equals the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo equals the third lowest double of 2/(transpose(v)*v);
 *   betahilolo equals the second lowest double of 2/(transpose(v)*v);
 *   betalololo equals the lowest double of 2/(transpose(v)*v). */

void CPU_cmplx8_factors_house 
 ( int n,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
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
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n         dimension of the vector x;
 *   xrehihihi has the highest doubles of the real parts of x;
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
 *   ximlololo has the lowest doubles of the imaginary parts of x;
 *   vrehihihi has space for n doubles;
 *   vrelohihi has space for n doubles;
 *   vrehilohi has space for n doubles;
 *   vrelolohi has space for n doubles;
 *   vrehihilo has space for n doubles;
 *   vrelohilo has space for n doubles;
 *   vrehilolo has space for n doubles;
 *   vrelololo has space for n doubles;
 *   vimhihihi has space for n doubles;
 *   vimlohihi has space for n doubles;
 *   vimhilohi has space for n doubles;
 *   vimlolohi has space for n doubles;
 *   vimhihilo has space for n doubles;
 *   vimlohilo has space for n doubles;
 *   vimhilolo has space for n doubles;
 *   vimlololo has space for n doubles.
 *
 * ON RETURN :
 *   vrehihihi has the highest doubles of the real parts
 *             of the Householder vector v;
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
 *   betalolohi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v). */

void CPU_dbl8_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double betahihihi, double betalohihi,
   double betahilohi, double betalolohi,
   double betahihilo, double betalohilo,
   double betahilolo, double betalololo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   Rhihihi  are the highest doubles of R;
 *   Rlohihi  are the second highest doubles of R;
 *   Rhilohi  are the third highest doubles of R;
 *   Rlolohi  are the fourth highest doubles of R;
 *   Rhihilo  are the fourth lowest doubles of R;
 *   Rlohilo  are the third lowest doubles of R;
 *   Rhilolo  are the second lowest doubles of R;
 *   Rlololo  are the lowest doubles of R;
 *   vhihihi  are the highest doubles of the Householder vector v;
 *   vlohihi  are the second highest doubles of v; 
 *   vhilohi  are the third highest doubles of v; 
 *   vlolohi  are the fourth highest doubles of v; 
 *   vhihilo  are the fourth lowest doubles of v;
 *   vlohilo  are the third lowest doubles of v;
 *   vhilolo  are the second lowest doubles of v;
 *   vlololo  are the lowest doubles of v;
 *   betahihihi is the highest double of 2/(transpose(v)*v);
 *   betalohihi is the second highest double of 2/(transpose(v)*v);
 *   betahilohi is the third highest double of 2/(transpose(v)*v);
 *   betalolohi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v).
 *
 * ON RETURN :
 *   Rhihihi  are the highest doubles of the update
 *            with the Householder matrix;
 *   Rlohihi  are the second highest doubles of the update;
 *   Rhilohi  are the third highest doubles of the update;
 *   Rlolohi  are the fourth highest doubles of the update;
 *   Rhihilo  are the fourth lowest doubles of the update;
 *   Rlohilo  are the third lowest doubles of the update;
 *   Rhilolo  are the second lowest doubles of the update;
 *   Rlololo  are the lowest doubles of the update. */

void CPU_cmplx8_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double betahihihi, double betalohihi,
   double betahilohi, double betalolohi,
   double betahihilo, double betalohilo,
   double betahilolo, double betalololo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows     number of rows in the matrix R;
 *   ncols     number of columns in the matrix R;
 *   k         current column index in R;
 *   Rrehihihi are the highest doubles of the real parts of R;
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
 *   vrehihihi are the highest doubles of the real parts 
 *             of the Householder vector v;
 *   vrelohihi are the second highest doubles of the real parts of v;
 *   vrehilohi are the third highest doubles of the real parts of v;
 *   vrelolohi are the fourth highest doubles of the real parts of v;
 *   vrehihilo are the fourth lowest doubles of the real parts of v;
 *   vrelohilo are the third lowest doubles of the real parts of v;
 *   vrehilolo are the second lowest doubles of the real parts of v;
 *   vrelololo are the lowest doubles of the real parts of v;
 *   vimhihihi are the highest doubles of the imaginary parts of v;
 *   vimlohihi are the second highest doubles of the imaginary parts of v;
 *   vimhilohi are the third highest doubles of the imaginary parts of v;
 *   vimlolohi are the fourth highest doubles of the imaginary parts of v;
 *   vimhihilo are the fourth lowest doubles of the imaginary parts of v;
 *   vimlohilo are the third lowest doubles of the imaginary parts of v;
 *   vimhilolo are the second lowest doubles of the imaginary parts of v;
 *   vimlololo are the lowest doubles of the imaginary parts of v;
 *   betahihihi is the highest double of 2/(transpose(v)*v);
 *   betalohihi is the second highest double of 2/(transpose(v)*v);
 *   betahilohi is the third highest double of 2/(transpose(v)*v);
 *   betalolohi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v).
 *
 * ON RETURN :
 *   Rrehihihi are the highest doubles of the real parts of the update
 *             with the Householder matrix;
 *   Rrelohihi are the second highest doubles of the real parts of the update;
 *   Rrehilohi are the third highest doubles of the real parts of the update;
 *   Rrelolohi are the fourth highest doubles of the real parts of the update;
 *   Rrehihilo are the fourth lowest doubles of the real parts of the update;
 *   Rrelohilo are the third lowest doubles of the real parts of the update;
 *   Rrehilolo are the second lowest doubles of the real parts of the update;
 *   Rrelololo are the lowest doubles of the real parts of the update;
 *   Rimhihihi are the highest doubles of the imaginary parts of the update;
 *   Rimlohihi are the second highest doubles of the imaginary parts
 *             of the update;
 *   Rimhilohi are the third highest doubles of the imaginary parts
 *             of the update;
 *   Rimlolohi are the fourth highest doubles of the imaginary parts
 *             of the update;
 *   Rimhihilo are the fourth lowest doubles of the imaginary parts
 *             of the update;
 *   Rimlohilo are the third lowest doubles of the imaginary parts
 *             of the update;
 *   Rimhilolo are the second lowest doubles of the imaginary parts
 *             of the update;
 *   Rimlololo are the lowest doubles of the imaginary parts of the update. */

void CPU_dbl8_factors_rightQupdate
 ( int n, int k,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double betahihihi, double betalohihi,
   double betahilohi, double betalolohi,
   double betahihilo, double betalohilo,
   double betahilolo, double betalololo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Qhihihi  are the highest doubles of Q;
 *   Qlohihi  are the second highest doubles of Q;
 *   Qhilohi  are the third highest doubles of Q;
 *   Qlolohi  are the fourth highest doubles of Q;
 *   Qhihilo  are the fourth lowest doubles of Q;
 *   Qlohilo  are the third lowest doubles of Q;
 *   Qhilolo  are the second lowest doubles of Q;
 *   Qlololo  are the lowest doubles of Q;
 *   vhihihi  are the highest doubles of the Householder vector v;
 *   vlohihi  are the second highest doubles of v;
 *   vhilohi  are the third highest doubles of v;
 *   vlolohi  are the fourth highest doubles of v;
 *   vhihilo  are the fourth lowest doubles of v;
 *   vlohilo  are the third lowest doubles of v;
 *   vhilolo  are the second lowest doubles of v;
 *   vlololo  are the lowest doubles of v;
 *   betahihihi is the highest double of 2/(transpose(v)*v);
 *   betalohihi is the second highest double of 2/(transpose(v)*v);
 *   betahilohi is the third highest double of 2/(transpose(v)*v);
 *   betalolohi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v).
 *
 * ON RETURN :
 *   Qhihihi  are the highest doubles of the update 
 *            with the Householder matrix;
 *   Qlohihi  are the second highest doubles of the update;
 *   Qhilohi  are the third highest doubles of the update;
 *   Qlolohi  are the fourth highest doubles of the update;
 *   Qhihilo  are the fourth lowest doubles of the update;
 *   Qlohilo  are the third lowest doubles of the update;
 *   Qhilolo  are the second lowest doubles of the update;
 *   Qlololo  are the lowest doubles of the update. */

void CPU_cmplx8_factors_rightQupdate
 ( int n, int k,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double betahihihi, double betalohihi,
   double betahilohi, double betalolohi,
   double betahihilo, double betalohilo,
   double betahilolo, double betalololo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n         dimension of the matrix Q;
 *   k         current column index in Q;
 *   Qrehihihi are the highest doubles of the real parts of Q;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third lowest doubles of the real parts of Q;
 *   Qrelolohi are the fourth lowest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo are the third lowest doubles of the real parts of Q;
 *   Qrehilolo are the second lowest doubles of the real parts of Q;
 *   Qrelololo are the lowest doubles of the real parts of Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the second highest doubles of the imaginary parts of Q;
 *   Qimlolohi are the second highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the second lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q;
 *   vrehihihi are the highest doubles of the real parts
 *             of the Householder vector v;
 *   vrelohihi are the second highest doubles of the real parts of v;
 *   vrehilohi are the third highest doubles of the real parts of v;
 *   vrelolohi are the fourth highest doubles of the real parts of v;
 *   vrehihilo are the fourth lowest doubles of the real parts of v;
 *   vrelohilo are the third lowest doubles of the real parts of v;
 *   vrehilolo are the second lowest doubles of the real parts of v;
 *   vrelololo are the lowest doubles of the real parts of v;
 *   vimhihihi are the highest doubles of the imaginary parts of v;
 *   vimlohihi are the second highest doubles of the imaginary parts of v;
 *   vimhilohi are the third highest doubles of the imaginary parts of v;
 *   vimlolohi are the fourth highest doubles of the imaginary parts of v;
 *   vimhihilo are the fourth lowest doubles of the imaginary parts of v;
 *   vimlohilo are the third lowest doubles of the imaginary parts of v;
 *   vimhilolo are the second lowest doubles of the imaginary parts of v;
 *   vimlololo are the lowest doubles of the imaginary parts of v;
 *   betahihihi is the highest double of 2/(transpose(v)*v);
 *   betalohihi is the second highest double of 2/(transpose(v)*v);
 *   betahilohi is the third highest double of 2/(transpose(v)*v);
 *   betalolohi is the fourth highest double of 2/(transpose(v)*v);
 *   betahihilo is the fourth lowest double of 2/(transpose(v)*v);
 *   betalohilo is the third lowest double of 2/(transpose(v)*v);
 *   betahilolo is the second lowest double of 2/(transpose(v)*v);
 *   betalololo is the lowest double of 2/(transpose(v)*v).
 *
 * ON RETURN :
 *   Qrehihihi has the highest doubles of the real parts of the update
 *             with the Householder matrix;
 *   Qrelohihi has the second highest doubles of the real parts of the update;
 *   Qrehilohi has the third highest doubles of the real parts of the update;
 *   Qrelolohi has the fourth highest doubles of the real parts of the update;
 *   Qrehihilo has the fourth lowest doubles of the real parts of the update;
 *   Qrelohilo has the third lowest doubles of the real parts of the update;
 *   Qrehilolo has the second lowest doubles of the real parts of the update;
 *   Qrelololo has the lowest doubles of the real parts of the update;
 *   Qimhihihi has the highest doubles of the imaginary parts of the update;
 *   Qimlohihi has the second highest doubles of the imaginary parts
 *             of the update;
 *   Qimhilohi has the third highest doubles of the imaginary parts
 *             of the update;
 *   Qimlolohi has the fourth highest doubles of the imaginary parts
 *             of the update;
 *   Qimhihilo has the fourth lowest doubles of the imaginary parts
 *             of the update;
 *   Qimlohilo has the third lowest doubles of the imaginary parts
 *             of the update;
 *   Qimhilolo has the second lowest doubles of the imaginary parts
 *             of the update;
 *   Qimlololo has the lowest doubles of the imaginary parts of the update. */

void CPU_dbl8_factors_houseqr
 ( int nrows, int ncols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo );
/*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   Ahihihi  has the highest doubles of A;
 *   Alohihi  has the second highest doubles of A;
 *   Ahilohi  has the third highest doubles of A;
 *   Alolohi  has the fourth highest doubles of A;
 *   Ahihilo  has the fourth lowest doubles of A;
 *   Alohilo  has the third lowest doubles of A;
 *   Ahilolo  has the second lowest doubles of A;
 *   Alololo  has the lowest doubles of A;
 *   Qhihihi  has space for an nrows-by-nrows matrix;
 *   Qlohihi  has space for an nrows-by-nrows matrix;
 *   Qhilohi  has space for an nrows-by-nrows matrix;
 *   Qlolohi  has space for an nrows-by-nrows matrix;
 *   Qhihilo  has space for an nrows-by-nrows matrix;
 *   Qlohilo  has space for an nrows-by-nrows matrix;
 *   Qhilolo  has space for an nrows-by-nrows matrix;
 *   Qlololo  has space for an nrows-by-nrows matrix;
 *   Rhihihi  has space for an nrows-by-ncols matrix;
 *   Rlohihi  has space for an nrows-by-ncols matrix;
 *   Rhilohi  has space for an nrows-by-ncols matrix;
 *   Rlolohi  has space for an nrows-by-ncols matrix.
 *   Rhihilo  has space for an nrows-by-ncols matrix;
 *   Rlohilo  has space for an nrows-by-ncols matrix;
 *   Rhilolo  has space for an nrows-by-ncols matrix;
 *   Rlololo  has space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qhihihi  are the highest doubles of an orthogonal matrix,
 *            transpose(Q)*A = R;
 *   Qlohihi  are the second highest doubles of Q
 *   Qhilohi  are the third highest doubles of Q
 *   Qlolohi  are the fourth highest doubles of Q
 *   Qhihilo  are the fourth lowest doubles of Q;
 *   Qlohilo  are the third lowest doubles of Q;
 *   Qhilolo  are the second lowest doubles of Q;
 *   Qlololo  are the lowest doubles of Q;
 *   Rhihihi  are the highest doubles of R;
 *   Rlohihi  are the second highest doubles of R;
 *   Rhilohi  are the third highest doubles of R;
 *   Rlolohi  are the fourth highest doubles of R;
 *   Rhihilo  are the fourth lowest doubles of R;
 *   Rlohilo  are the third lowest doubles of R;
 *   Rhilolo  are the second lowest doubles of R;
 *   Rlololo  are the lowest doubles of R. */

void CPU_cmplx8_factors_houseqr
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
   double **Rimhilolo, double **Rimlololo );
/*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows     number of rows of A;
 *   ncols     number of columns of A;
 *   Arehihihi are the highest doubles of the real parts of A;
 *   Arelohihi are the second highest doubles of the real parts of A;
 *   Arehilohi are the third highest doubles of the real parts of A;
 *   Arelolohi are the fourth highest doubles of the real parts of A;
 *   Arehihilo are the fourth lowest doubles of the real parts of :
 *   Arelohilo are the third lowest doubles of the real parts of :
 *   Arehilolo are the second lowest doubles of the real parts of :
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
 *   Rimhihihi has space for an nrows-by-ncols matrix;
 *   Rimlohihi has space for an nrows-by-ncols matrix;
 *   Rimhilolo has space for an nrows-by-ncols matrix;
 *   Rimlololo has space for an nrows-by-ncols matrix.
 *   Rimhilolo has space for an nrows-by-ncols matrix;
 *   Rimlololo has space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qrehihihi are the highest doubles of the real parts
 *             of the orthogonal matrix Q, transpose(Q)*A = R;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third highest doubles of the real parts of Q;
 *   Qrelolohi are the fourth highest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo are the third lowest doubles of the real parts of Q;
 *   Qrehilolo are the second lowest doubles of the real parts of Q;
 *   Qrelololo are the lowest doubles of the real parts of Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the third highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q;
 *   Rrehihihi are the highest doubles of the real parts of R,
 *             the reduced upper triangular form of A;
 *   Rrelohihi are the second highest doubles of the real parts of R,
 *   Rrehilohi are the third highest doubles of the real parts of R,
 *   Rrelolohi are the fourth highest doubles of the real parts of R,
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

void CPU_dbl8_factors_qrbs
 ( int nrows, int ncols,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *rhshihihi, double *rhslohihi, double *rhshilohi, double *rhslolohi,
   double *rhshihilo, double *rhslohilo, double *rhshilolo, double *rhslololo,
   double *solhihihi, double *sollohihi, double *solhilohi, double *sollolohi,
   double *solhihilo, double *sollohilo, double *solhilolo, double *sollololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo ); 
/*
 * DESCRIPTION :
 *   Applies back substitution for the least squares solution
 *   with the QR factorization.
 *
 * REQUIRED : nrows >= ncols;
 *
 * ON ENTRY :
 *   nrows    number of rows of R, dimension of Q;
 *   ncols    number of columns of R;
 *   Qhihihi  highest doubles of the Q of the QR factorization;
 *   Qlohihi  second highest doubles of the Q of the QR factorization;
 *   Qhilohi  third highest doubles of the Q of the QR factorization;
 *   Qlolohi  fourth highest doubles of the Q of the QR factorization;
 *   Qhihilo  fourth lowest doubles of the Q of the QR factorization;
 *   Qlohilo  third lowest doubles of the Q of the QR factorization;
 *   Qhilolo  second lowest doubles of the Q of the QR factorization;
 *   Qlololo  lowest doubles of the Q of the QR factorization;
 *   Rhihihi  highest doubles of the R of the QR factorization;
 *   Rlohihi  second highest doubles of the R of the QR factorization;
 *   Rhilohi  third highest doubles of the R of the QR factorization;
 *   Rlolohi  fourth highest doubles of the R of the QR factorization;
 *   Rhihilo  fourth lowest doubles of the R of the QR factorization;
 *   Rlohilo  third lowest doubles of the R of the QR factorization;
 *   Rhilolo  second lowest doubles of the R of the QR factorization;
 *   Rlololo  lowest doubles of the R of the QR factorization;
 *   rhshihihi has the highest doubles of the right hand side, of size nrows;
 *   rhslohihi has the second highest doubles of the right hand side;
 *   rhshilohi has the third highest doubles of the right hand side;
 *   rhslolohi has the fourth highest doubles of the right hand side;
 *   rhshihilo has the fourth lowest doubles of the right hand side;
 *   rhslohilo has the third lowest doubles of the right hand side;
 *   rhshilolo has the second lowest doubles of the right hand side;
 *   rhslololo has the lowest doubles of the right hand side;
 *   solhihihi has space for the highest doubles of the solution, of size ncols;
 *   sollohihi has space for the second highest doubles of the solution;
 *   solhilohi has space for the third highest doubles of the solution;
 *   sollolohi has space for the fourth highest doubles of the solution;
 *   solhihilo has space for the fourth lowest doubles of the solution;
 *   sollohilo has space for the third lowest doubles of the solution;
 *   solhilolo has space for the second lowest doubles of the solution;
 *   sollololo has space for the lowest doubles of the solution;
 *   wrkvechihihi has work space for nrows elements;
 *   wrkveclohihi has work space for nrows elements;
 *   wrkvechilohi has work space for nrows elements;
 *   wrkveclolohi has work space for nrows elements;
 *   wrkvechihilo has work space for nrows elements;
 *   wrkveclohilo has work space for nrows elements;
 *   wrkvechilolo has work space for nrows elements;
 *   wrkveclololo has work space for nrows elements.
 *
 * ON RETURN :
 *   solhihihi has the highest doubles of the least squares solution;
 *   sollohihi has the second highest doubles of the least squares solution;
 *   solhilohi has the third highest doubles of the least squares solution;
 *   sollolohi has the fourth highest doubles of the least squares solution;
 *   solhihilo has the fourth lowest doubles of the least squares solution;
 *   sollohilo has the third lowest doubles of the least squares solution;
 *   solhilolo has the second lowest doubles of the least squares solution;
 *   sollololo has the lowest doubles of the least squares solution;
 *   wrkvechihihi has the highest doubles of the product of Q^T with rhs;
 *   wrkveclohihi has the second highest doubles of Q^T*rhs;
 *   wrkvechilohi has the third highest doubles of Q^T*rhs;
 *   wrkveclolohi has the fourth highest doubles of Q^T*rhs;
 *   wrkvechihilo has the fourth lowest doubles of Q^T*rhs;
 *   wrkveclohilo has the third lowest doubles of Q^T*rhs;
 *   wrkvechilolo has the second lowest doubles of Q^T*rhs;
 *   wrkveclololo has the lowest doubles of Q^T*rhs. */

void CPU_cmplx8_factors_qrbs
 ( int nrows, int ncols,
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
   double *rhsrehihihi, double *rhsrelohihi,
   double *rhsrehilohi, double *rhsrelolohi,
   double *rhsrehihilo, double *rhsrelohilo,
   double *rhsrehilolo, double *rhsrelololo,
   double *rhsimhihihi, double *rhsimlohihi,
   double *rhsimhilohi, double *rhsimlolohi,
   double *rhsimhihilo, double *rhsimlohilo,
   double *rhsimhilolo, double *rhsimlololo,
   double *solrehihihi, double *solrelohihi,
   double *solrehilohi, double *solrelolohi,
   double *solrehihilo, double *solrelohilo,
   double *solrehilolo, double *solrelololo,
   double *solimhihihi, double *solimlohihi,
   double *solimhilohi, double *solimlolohi,
   double *solimhihilo, double *solimlohilo,
   double *solimhilolo, double *solimlololo,
   double *wrkvecrehihihi, double *wrkvecrelohihi,
   double *wrkvecrehilohi, double *wrkvecrelolohi,
   double *wrkvecrehihilo, double *wrkvecrelohilo,
   double *wrkvecrehilolo, double *wrkvecrelololo,
   double *wrkvecimhihihi, double *wrkvecimlohihi,
   double *wrkvecimhilohi, double *wrkvecimlolohi,
   double *wrkvecimhihilo, double *wrkvecimlohilo,
   double *wrkvecimhilolo, double *wrkvecimlololo );
/*
 * DESCRIPTION :
 *   Applies back substitution for the least squares solution
 *   with the QR factorization, on complex data.
 *
 * REQUIRED : nrows >= ncols;
 *
 * ON ENTRY :
 *   nrows    number of rows of R, dimension of Q;
 *   ncols    number of columns of R;
 *   Qrehihihi are the highest doubles of the real parts of Q of the QR;
 *   Qrelohihi are the 2nd highest doubles of the real parts of Q of the QR;
 *   Qrehilohi are the 3rd highest doubles of the real parts of Q of the QR;
 *   Qrelolohi are the 4th highest doubles of the real parts of Q of the QR;
 *   Qrehihilo are the 4th lowest doubles of the real parts of Q of the QR;
 *   Qrelohilo are the 3rd lowest doubles of the real parts of Q of the QR;
 *   Qrehilolo are the 2nd lowest doubles of the real parts of Q of the QR;
 *   Qrelololo are the lowest doubles of the real parts of Q of the QR;
 *   Qimhihihi are the highest doubles of the imag parts of Q of the QR;
 *   Qimlohihi are the 2nd highest doubles of the imag parts of Q of the QR;
 *   Qimhilohi are the 3rd highest doubles of the imag parts of Q of the QR;
 *   Qimlolohi are the 4th highest doubles of the imag parts of Q of the QR;
 *   Qimhihilo are the 4th lowest doubles of the imag parts of Q of the QR;
 *   Qimlohilo are the 3rd lowest doubles of the imag parts of Q of the QR;
 *   Qimhilolo are the 2nd lowest doubles of the imag parts of Q of the QR;
 *   Qimlololo are the lowest doubles of the imag parts of Q of the QR;
 *   Rrehihihi are the highest doubles of the real parts of R of the QR;
 *   Rrelohihi are the 2nd highest doubles of the real parts of R of the QR;
 *   Rrehilohi are the 3rd highest doubles of the real parts of R of the QR;
 *   Rrelolohi are the 4th highest doubles of the real parts of R of the QR;
 *   Rrehihilo are the 4th lowest doubles of the real parts of R of the QR;
 *   Rrelohilo are the 3rd lowest doubles of the real parts of R of the QR;
 *   Rrehilolo are the 2nd lowest doubles of the real parts of R of the QR;
 *   Rrelololo are the lowest doubles of the real parts of R of the QR;
 *   Rimhihihi are the highest doubles of the imag parts of R of the QR;
 *   Rimlohihi are the 2nd highest doubles of the imag parts of R of the QR;
 *   Rimhilohi are the 3rd highest doubles of the imag parts of R of the QR;
 *   Rimlolohi are the 4th highest doubles of the imag parts of R of the QR;
 *   Rimhihilo are the 4th lowest doubles of the imag parts of R of the QR;
 *   Rimlohilo are the 3rd lowest doubles of the imag parts of R of the QR;
 *   Rimhilolo are the 2nd lowest doubles of the imag parts of R of the QR;
 *   Rimlololo are the lowest doubles of the imag parts of R of the QR;
 *   rhsrehihihi are the highest doubles of the real parts of right hand side,
 *            rhs has size nrows;
 *   rhsrelohihi are the 2nd highest doubles of the real parts of rhs;
 *   rhsrehilohi are the 3rd highest doubles of the real parts of rhs;
 *   rhsrelolohi are the 4th highest doubles of the real parts of rhs;
 *   rhsrehihilo are the 4th lowest doubles of the real parts of rhs;
 *   rhsrelohilo are the 3rd lowest doubles of the real parts of rhs;
 *   rhsrehilolo are the 2nd lowest doubles of the real parts of rhs;
 *   rhsrelololo are the lowest doubles of the real parts of rhs;
 *   rhsimhihihi are the highest doubles of the imag parts of rhs;
 *   rhsimlohihi are the 2nd highest doubles of the imag parts of rhs;
 *   rhsimhilohi are the 3rd highest doubles of the imag parts of rhs;
 *   rhsimlolohi are the 4th highest doubles of the imag parts of rhs;
 *   rhsimhihilo are the 4th lowest doubles of the imag parts of rhs;
 *   rhsimlohilo are the 3rd lowest doubles of the imag parts of rhs;
 *   rhsimhilolo are the 2nd lowest doubles of the imag parts of rhs;
 *   rhsimlololo are lowest doubles of the imag parts of rhs;
 *   solrehihihi has space for the highest doubles of the real parts
 *            of the solution, of size ncols;
 *   solrelohihi has space for the second highest doubles of the real parts
 *            of the solution, of size ncols;
 *   solrehilohi has space for the third highest doubles of the real parts
 *            of the solution, of size ncols;
 *   solrelolohi has space for the fourth highest doubles of the real parts
 *            of the solution, of size ncols;
 *   solrehihilo has space for the fourth lowest doubles of the real parts
 *            of the solution;
 *   solrelohilo has space for the third lowest doubles of the real parts
 *            of the solution;
 *   solrehilolo has space for the second lowest doubles of the real parts
 *            of the solution;
 *   solrelololo has space for the lowest doubles of the real parts
 *            of the solution;
 *   solimhihihi has space for the highest doubles of the imaginary parts
 *            of the solution;
 *   solimlohihi has space for the second highest doubles of the imag parts
 *            of the solution;
 *   solimhilohi has space for the third highest doubles of the imag parts
 *            of the solution;
 *   solimlolohi has space for the fourth highest doubles of the imag parts
 *            of the solution;
 *   solimhihilo has space for the fourth lowest doubles of the imag parts
 *            of the solution;
 *   solimlohilo has space for the third lowest doubles of the imag parts
 *            of the solution;
 *   solimhilolo has space for the second lowest doubles of the imag parts
 *            of the solution;
 *   solimlololo has space for the lowest doubles of the imaginary parts
 *            of the solution;
 *   wrkvecrehihihi is work space for nrows elements;
 *   wrkvecrelohihi is work space for nrows elements;
 *   wrkvecrehilohi is work space for nrows elements;
 *   wrkvecrelolohi is work space for nrows elements;
 *   wrkvecrehihilo is work space for nrows elements;
 *   wrkvecrelohilo is work space for nrows elements;
 *   wrkvecrehilolo is work space for nrows elements;
 *   wrkvecrelololo is work space for nrows elements;
 *   wrkvecimhihihi is work space for nrows elements;
 *   wrkvecimlohihi is work space for nrows elements;
 *   wrkvecimhilohi is work space for nrows elements;
 *   wrkvecimlolohi is work space for nrows elements.
 *   wrkvecimhihilo is work space for nrows elements;
 *   wrkvecimlohilo is work space for nrows elements;
 *   wrkvecimhilolo is work space for nrows elements;
 *   wrkvecimlololo is work space for nrows elements.
 *
 * ON RETURN :
 *   solrehihihi are the highest doubles of the real parts of the solution;
 *   solrelohihi are the 2nd highest doubles of the real parts of the solution;
 *   solrehilohi are the 2rd highest doubles of the real parts of the solution;
 *   solrelolohi are the 4th highest doubles of the real parts of the solution;
 *   solrehihilo are the 4th lowest doubles of the real parts of the solution;
 *   solrelohilo are the 3rd lowest doubles of the real parts of the solution;
 *   solrehilolo are the 2nd lowest doubles of the real parts of the solution;
 *   solrelololo are the lowest doubles of the real parts of the solution;
 *   solimhihihi are the highest doubles of the imag parts of the solution;
 *   solimlohihi are the 2nd highest doubles of the imag parts of sol;
 *   solimhilohi are the 3rd highest doubles of the imag parts of sol;
 *   solimlolohi are the 4th highest doubles of the imag parts of sol;
 *   solimhihilo are the 4th lowest doubles of the imag parts of sol;
 *   solimlohilo are the 3rd lowest doubles of the imag parts of sol;
 *   solimhilolo are the 2nd lowest doubles of the imag parts of sol;
 *   solimlololo are the lowest doubles of the imag parts of the solution;
 *   wrkvecrehihihi has the highest doubles of the real parts of Q^H*rhs;
 *   wrkvecrelohihi has the 2nd highest doubles of the real parts of Q^H*rhs;
 *   wrkvecrehilohi has the 3rd highest doubles of the real parts of Q^H*rhs;
 *   wrkvecrelolohi has the 4th highest doubles of the real parts of Q^H*rhs;
 *   wrkvecrehihilo has the 4th lowest doubles of the real parts of Q^H*rhs;
 *   wrkvecrelohilo has the 3rd lowest doubles of the real parts of Q^H*rhs;
 *   wrkvecrehilolo has the 2nd lowest doubles of the real parts of Q^H*rhs;
 *   wrkvecrelololo has the lowest doubles of the real parts of Q^H*rhs;
 *   wrkvecimhihihi has the highest doubles of the imaginary parts of Q^H*rhs;
 *   wrkvecimlohihi has the second highest doubles of the imaginary parts
 *            of Q^H*rhs;
 *   wrkvecimhilohi has the third highest doubles of the imaginary parts
 *            of Q^H*rhs;
 *   wrkvecimlolohi has the fourth highest doubles of the imaginary parts
 *            of Q^H*rhs;
 *   wrkvecimhihilo has the fourth lowest doubles of the imaginary parts 
 *            of Q^H*rhs;
 *   wrkvecimlohilo has the third lowest doubles of the imaginary parts 
 *            of Q^H*rhs;
 *   wrkvecimhilolo has the second lowest doubles of the imaginary parts 
 *            of Q^H*rhs;
 *   wrkvecimlololo has the lowest doubles of the imaginary parts
 *            of Q^H*rhs. */

#endif
