with Standard_Floating_Ring;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Abstract_Ring.Field;

package Standard_Floating_Ring.FField is
  new Standard_Floating_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of standard floating-point numbers
--   to a field.
