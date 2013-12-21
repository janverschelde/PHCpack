with Multprec_Floating_Ring;
with Multprec_Floating_Ring.FField;
with Generic_Complex_Numbers;

package Multprec_Complex_Numbers is 
  new Generic_Complex_Numbers(Multprec_Floating_Ring,
                              Multprec_Floating_Ring.FField);

-- DESCRIPTION :
--   Defines the multi-precision complex numbers.
