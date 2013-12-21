with Multprec_Floating_Ring;
with Multprec_Floating_Ring.Field;
with Multprec_Complex_Ring;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Abstract_Ring.Complex;

package Multprec_Complex_Ring.Complex is
  new Multprec_Complex_Ring.Complex(Multprec_Floating_Ring,
                                    Multprec_Floating_Ring.Field,
                                    AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of multi-precision complex numbers
--   to a complex field.
