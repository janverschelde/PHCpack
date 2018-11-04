with DoblDobl_Complex_Series_io;
with DoblDobl_Complex_Series_Ring;
with Abstract_Ring_io;

package DoblDobl_Complex_Series_Ring_io is
  new Abstract_Ring_io(DoblDobl_Complex_Series_Ring,
                       DoblDobl_Complex_Series_io.get,
                       DoblDobl_Complex_Series_io.put,
                       DoblDobl_Complex_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in double double precision.
