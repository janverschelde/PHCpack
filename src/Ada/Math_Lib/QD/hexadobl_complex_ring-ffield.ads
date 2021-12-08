with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Numbers;           use HexaDobl_Complex_Numbers;
with Abstract_Ring.Field;

package HexaDobl_Complex_Ring.FField is
  new HexaDobl_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of complex hexa double numbers to a field.
