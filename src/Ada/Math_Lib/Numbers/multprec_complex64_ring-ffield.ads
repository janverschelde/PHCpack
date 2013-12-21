with Multprec_Complex64_Ring;
with Multprec_Complex64_Numbers;         use Multprec_Complex64_Numbers;
with Abstract_Ring.Field;

package Multprec_Complex64_Ring.FField is
  new Multprec_Complex64_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of multi-precision complex numbers
--   to a field.
