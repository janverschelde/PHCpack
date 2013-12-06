with Standard_Complex_Ring;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Abstract_Ring.Field;

package Standard_Complex_Ring.FField is
  new Standard_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of standard complex numbers to a field.
