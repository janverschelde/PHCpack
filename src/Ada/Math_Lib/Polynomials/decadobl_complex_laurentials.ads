with Generic_Laurent_Polynomials;
with DecaDobl_Complex_Ring;

package DecaDobl_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(DecaDobl_Complex_Ring);

-- DESCRIPTION :
--   Defines Laurent polynomials over the ring of complex deca doubles.
--   The "Laurential" is a contraction of "Laurent polynomial".
