with Generic_Laurent_Polynomials;
with OctoDobl_Complex_Ring;

package OctoDobl_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(OctoDobl_Complex_Ring);

-- DESCRIPTION :
--   Defines Laurent polynomials over the ring of complex octo doubles.
--   The "Laurential" is a contraction of "Laurent polynomial".
