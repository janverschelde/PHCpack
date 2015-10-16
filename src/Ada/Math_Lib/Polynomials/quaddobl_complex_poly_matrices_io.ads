with QuadDobl_Complex_Poly_Ring_io;
with QuadDobl_Complex_Poly_Vectors;
with QuadDobl_Complex_Poly_Matrices;
with Generic_Matrices_io;

package QuadDobl_Complex_Poly_Matrices_io is 
  new Generic_Matrices_io(QuadDobl_Complex_Poly_Ring_io,
                          QuadDobl_Complex_Poly_Vectors,
                          QuadDobl_Complex_Poly_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of complex polynomials,
--   with coefficients in quad double precision.
