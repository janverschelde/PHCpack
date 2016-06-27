with QuadDobl_Dense_Series_Ring;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Series_Polynomials;
with Generic_Polynomial_Functions;

package QuadDobl_Series_Poly_Functions is
  new Generic_Polynomial_Functions(QuadDobl_Dense_Series_Ring,
                                   QuadDobl_Dense_Series_Vectors,
                                   QuadDobl_Series_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with quad double complex coefficients.
