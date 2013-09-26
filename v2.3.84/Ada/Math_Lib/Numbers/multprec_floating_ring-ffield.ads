with Multprec_Floating_Ring;
with Multprec_Floating_Numbers;           use Multprec_Floating_Numbers;
with Abstract_Ring.Field;

package Multprec_Floating_Ring.FField is
  new Multprec_Floating_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of multi-precision floating-point
--   numbers to a field.
