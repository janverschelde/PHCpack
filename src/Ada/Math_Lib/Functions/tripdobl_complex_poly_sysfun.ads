with TripDobl_Complex_Ring;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Functions;
with TripDobl_Complex_Poly_Systems;
with Generic_Poly_System_Functions;

package TripDobl_Complex_Poly_SysFun is
  new Generic_Poly_System_Functions(TripDobl_Complex_Ring,
                                    TripDobl_Complex_Vectors,
                                    TripDobl_Complex_VecVecs,
                                    TripDobl_Complex_Polynomials,
                                    TripDobl_Complex_Poly_Functions,
                                    TripDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials over the
--   ring of triple double complex numbers.
