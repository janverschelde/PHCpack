with DecaDobl_Complex_Series_Ring;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_CSeries_Polynomials;
with Generic_Polynomial_Functions;

package DecaDobl_CSeries_Poly_Functions is
  new Generic_Polynomial_Functions(DecaDobl_Complex_Series_Ring,
                                   DecaDobl_Complex_Series_Vectors,
                                   DecaDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with deca double complex coefficients.
