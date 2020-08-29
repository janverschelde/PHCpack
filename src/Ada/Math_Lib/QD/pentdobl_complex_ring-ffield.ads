with PentDobl_Complex_Ring;
with PentDobl_Complex_Numbers;           use PentDobl_Complex_Numbers;
with Abstract_Ring.Field;

package PentDobl_Complex_Ring.FField is
  new PentDobl_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of complex penta double numbers to a field.
