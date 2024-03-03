with TripDobl_Complex_Ring;
with TripDobl_Complex_Polynomials;
with Generic_Lists_of_Terms;

package TripDobl_Complex_Term_Lists is
  new Generic_Lists_of_Terms(TripDobl_Complex_Ring,
                             TripDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms for polynomial with complex coefficients 
--   in triple double precision.
