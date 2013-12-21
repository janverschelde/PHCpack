with Generic_Laurent_Polynomials;
with DoblDobl_Complex_Ring;

package DoblDobl_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(DoblDobl_Complex_Ring);

-- DESCRIPTION :
--   Defines Laurent polynomials over the ring of complex double doubles.
--   The "Laurential" is a contraction of "Laurent polynomial".
