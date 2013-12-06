with Multprec_Floating64_Ring;
with Multprec_Floating64_Numbers;         use Multprec_Floating64_Numbers;
with Abstract_Ring.Field;

package Multprec_Floating64_Ring.FField is
  new Multprec_Floating64_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of multi-precision floating-point
--   numbers to a field.
