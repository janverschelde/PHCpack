with TripDobl_Complex_Series_Ring;
with TripDobl_CSeries_Polynomials;
with Generic_Polynomial_Systems;

package TripDobl_CSeries_Poly_Systems is
  new Generic_Polynomial_Systems(TripDobl_Complex_Series_Ring,
                                 TripDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines systems of polynomials in several variables with coefficients
--   as series of triple double precision complex numbers.
