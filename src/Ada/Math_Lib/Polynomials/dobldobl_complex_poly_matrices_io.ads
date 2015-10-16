with DoblDobl_Complex_Poly_Ring_io;
with DoblDobl_Complex_Poly_Vectors;
with DoblDobl_Complex_Poly_Matrices;
with Generic_Matrices_io;

package DoblDobl_Complex_Poly_Matrices_io is 
  new Generic_Matrices_io(DoblDobl_Complex_Poly_Ring_io,
                          DoblDobl_Complex_Poly_Vectors,
                          DoblDobl_Complex_Poly_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of complex polynomials,
--   with coefficients in double double precision.
