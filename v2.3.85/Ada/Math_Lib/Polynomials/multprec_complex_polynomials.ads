with Generic_Polynomials;
with Multprec_Complex_Ring;

package Multprec_Complex_Polynomials is 
  new Generic_Polynomials(Multprec_Complex_Ring);

-- DESCRIPTION :
--   Defines the polynomials over the ring of multi-precision complex numbers.
