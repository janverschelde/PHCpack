with TripDobl_Complex_Series_Ring;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Series_VecVecs;
with TripDobl_Complex_Series_Matrices;
with TripDobl_CSeries_Polynomials;
with TripDobl_CSeries_Poly_Functions;
with TripDobl_CSeries_Poly_Systems;
with TripDobl_CSeries_Poly_SysFun;
with Generic_Jacobian_Matrices;

package TripDobl_CSeries_Jaco_Matrices is
  new Generic_Jacobian_Matrices(TripDobl_Complex_Series_Ring,
                                TripDobl_Complex_Series_Vectors,
                                TripDobl_Complex_Series_VecVecs,
                                TripDobl_Complex_Series_Matrices,
                                TripDobl_CSeries_Polynomials,
                                TripDobl_CSeries_Poly_Functions,
                                TripDobl_CSeries_Poly_Systems,
                                TripDobl_CSeries_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials in several variables over the field of
--   truncated power series with triple double precision complex numbers.
