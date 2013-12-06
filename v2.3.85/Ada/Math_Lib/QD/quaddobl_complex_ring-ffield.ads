with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Abstract_Ring.Field;

package QuadDobl_Complex_Ring.FField is
  new QuadDobl_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of quad double numbers
--   to a field.
