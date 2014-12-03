with Multprec_Complex_Ring;
with Multprec_Complex_Polynomials;
with Generic_Lists_of_Terms;

package Multprec_Complex_Term_Lists is
  new Generic_Lists_of_Terms(Multprec_Complex_Ring,
                             Multprec_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms for polynomial with complex coefficients 
--   in arbitrary multiprecision.
