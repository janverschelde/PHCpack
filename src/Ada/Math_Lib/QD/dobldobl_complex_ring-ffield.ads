with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Abstract_Ring.Field;

package DoblDobl_Complex_Ring.FField is
  new DoblDobl_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of double double numbers
--   to a field.
