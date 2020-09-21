with TripDobl_Complex_Series_Ring;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Series_VecVecs;
with TripDobl_CSeries_Polynomials;
with TripDobl_CSeries_Poly_Functions;
with TripDobl_CSeries_Poly_Systems;
with Generic_Poly_System_Functions;

package TripDobl_CSeries_Poly_SysFun is
  new Generic_Poly_System_Functions(TripDobl_Complex_Series_Ring,
                                    TripDobl_Complex_Series_Vectors,
                                    TripDobl_Complex_Series_VecVecs,
                                    TripDobl_CSeries_Polynomials,
                                    TripDobl_CSeries_Poly_Functions,
                                    TripDobl_CSeries_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   triple double precision complex numbers as coefficients.
