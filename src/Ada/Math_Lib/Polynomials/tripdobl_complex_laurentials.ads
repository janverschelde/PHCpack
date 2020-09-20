with Generic_Laurent_Polynomials;
with TripDobl_Complex_Ring;

package TripDobl_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(TripDobl_Complex_Ring);

-- DESCRIPTION :
--   Defines Laurent polynomials over the ring of complex triple doubles.
--   The "Laurential" is a contraction of "Laurent polynomial".
