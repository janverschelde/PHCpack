with TripDobl_Complex_Ring;
with TripDobl_Complex_Numbers;           use TripDobl_Complex_Numbers;
with Abstract_Ring.Field;

package TripDobl_Complex_Ring.FField is
  new TripDobl_Complex_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of complex triple double numbers to a field.
