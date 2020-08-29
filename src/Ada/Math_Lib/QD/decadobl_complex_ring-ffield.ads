with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Numbers;           use DecaDobl_Complex_Numbers;
with Abstract_Ring.Field;

package DecaDobl_Complex_Ring.FField is
  new DecaDobl_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of complex deca double numbers to a field.
