with Generic_Polynomials;
with Multprec_Floating_Ring;

package Multprec_Floating_Polynomials is 
  new Generic_Polynomials(Multprec_Floating_Ring);

-- DESCRIPTION :
--   Defines polynomials with multiprecision floating-point coefficients.
