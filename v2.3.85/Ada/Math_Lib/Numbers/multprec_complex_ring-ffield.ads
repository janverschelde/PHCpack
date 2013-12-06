with Multprec_Complex_Ring;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Abstract_Ring.Field;

package Multprec_Complex_Ring.FField is
  new Multprec_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of multi-precision complex numbers
--   to a field.
