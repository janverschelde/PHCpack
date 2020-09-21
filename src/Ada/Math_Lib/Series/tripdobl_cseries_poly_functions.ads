with TripDobl_Complex_Series_Ring;
with TripDobl_Complex_Series_Vectors;
with TripDobl_CSeries_Polynomials;
with Generic_Polynomial_Functions;

package TripDobl_CSeries_Poly_Functions is
  new Generic_Polynomial_Functions(TripDobl_Complex_Series_Ring,
                                   TripDobl_Complex_Series_Vectors,
                                   TripDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with triple double complex coefficients.
