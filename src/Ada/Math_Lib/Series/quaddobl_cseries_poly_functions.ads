with QuadDobl_Complex_Series_Ring;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Polynomials;
with Generic_Polynomial_Functions;

package QuadDobl_CSeries_Poly_Functions is
  new Generic_Polynomial_Functions(QuadDobl_Complex_Series_Ring,
                                   QuadDobl_Complex_Series_Vectors,
                                   QuadDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with quad double complex coefficients.
