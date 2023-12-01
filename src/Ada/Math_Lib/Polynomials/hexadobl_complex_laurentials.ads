with Generic_Laurent_Polynomials;
with HexaDobl_Complex_Ring;

package HexaDobl_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(HexaDobl_Complex_Ring);

-- DESCRIPTION :
--   Defines Laurent polynomials over the ring of complex hexa doubles.
--   The "Laurential" is a contraction of "Laurent polynomial".
