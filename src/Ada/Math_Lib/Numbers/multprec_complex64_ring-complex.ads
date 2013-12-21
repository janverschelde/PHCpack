with Multprec_Floating64_Ring;
with Multprec_Floating64_Ring.Field;
with Multprec_Complex64_Ring;
with Multprec_Complex64_Numbers;         use Multprec_Complex64_Numbers;
with Abstract_Ring.Complex;

package Multprec_Complex64_Ring.Complex is
  new Multprec_Complex64_Ring.Complex(Multprec_Floating64_Ring,
                                      Multprec_Floating64_Ring.Field,
                                      AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of multi-precision complex numbers
--   to a complex field.
