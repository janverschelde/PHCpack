with Generic_Laurent_Polynomials;
with Multprec_Complex_Ring;

package Multprec_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(Multprec_Complex_Ring);

-- DESCRIPTION :
--   Defines the Laurent polynomials with multi-precision complex
--   coefficients.  The "Laurential" is a contraction of Laurent polynomial.
