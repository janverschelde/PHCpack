with OctoDobl_Complex_Series_Ring;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_CSeries_Polynomials;
with Generic_Polynomial_Functions;

package OctoDobl_CSeries_Poly_Functions is
  new Generic_Polynomial_Functions(OctoDobl_Complex_Series_Ring,
                                   OctoDobl_Complex_Series_Vectors,
                                   OctoDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with octo double complex coefficients.
