with Generic_Laurent_Polynomials;
with PentDobl_Complex_Ring;

package PentDobl_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(PentDobl_Complex_Ring);

-- DESCRIPTION :
--   Defines Laurent polynomials over the ring of complex penta doubles.
--   The "Laurential" is a contraction of "Laurent polynomial".
