with Generic_Laurent_Polynomials;
with QuadDobl_Complex_Ring;

package QuadDobl_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(QuadDobl_Complex_Ring);

-- DESCRIPTION :
--   Defines Laurent polynomials over the ring of complex quad doubles.
--   The "Laurential" is a contraction of "Laurent polynomial".
