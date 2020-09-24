with TripDobl_Complex_Ring;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Polynomials;
with Generic_Polynomial_Functions;

package TripDobl_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(TripDobl_Complex_Ring,
                                   TripDobl_Complex_Vectors,
                                   TripDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for complex triple double numbers.
