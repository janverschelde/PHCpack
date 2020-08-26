with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Numbers;           use OctoDobl_Complex_Numbers;
with Abstract_Ring.Field;

package OctoDobl_Complex_Ring.FField is
  new OctoDobl_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of complex octo double numbers to a field.
