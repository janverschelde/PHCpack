with TripDobl_Complex_Ring;
with TripDobl_Complex_Polynomials;
with Generic_Polynomial_Systems;

package TripDobl_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(TripDobl_Complex_Ring,
                                 TripDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the triple double complex numbers.
