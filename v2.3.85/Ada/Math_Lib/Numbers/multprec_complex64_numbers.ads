with Multprec_Floating64_Ring;
with Multprec_Floating64_Ring.FField;
with Generic_Complex_Numbers;

package Multprec_Complex64_Numbers is 
  new Generic_Complex_Numbers(Multprec_Floating64_Ring,
                              Multprec_Floating64_Ring.FField);

-- DESCRIPTION :
--   Defines the multi-precision complex numbers over the field of
--   floating-point numbers with 64-bit fraction and exponent arithmetic.
