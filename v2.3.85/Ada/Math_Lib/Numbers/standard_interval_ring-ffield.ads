with Standard_Interval_Ring;
with Standard_Interval_Numbers;           use Standard_Interval_Numbers;
with Abstract_Ring.Field;

package Standard_Interval_Ring.FField is
  new Standard_Interval_Ring.Field("<",">",Width,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of standard interval numbers
--   to a field.
