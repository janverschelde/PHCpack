with Generic_Polynomials;
with TripDobl_Complex_Series_Ring;

package TripDobl_CSeries_Polynomials is 
  new Generic_Polynomials(TripDobl_Complex_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with triple double complex coefficients.
