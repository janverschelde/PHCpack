with TripDobl_Complex_Series_io;
with TripDobl_Complex_Series_Ring;
with Abstract_Ring_io;

package TripDobl_Complex_Series_Ring_io is
  new Abstract_Ring_io(TripDobl_Complex_Series_Ring,
                       TripDobl_Complex_Series_io.get,
                       TripDobl_Complex_Series_io.put,
                       TripDobl_Complex_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in triple double precision.
