with TripDobl_Complex_Numbers_io;
with TripDobl_Complex_Ring;
with Abstract_Ring_io;

package TripDobl_Complex_Ring_io is
  new Abstract_Ring_io(TripDobl_Complex_Ring,
                       TripDobl_Complex_Numbers_io.get,
                       TripDobl_Complex_Numbers_io.put,
                       TripDobl_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for triple double complex numbers.
