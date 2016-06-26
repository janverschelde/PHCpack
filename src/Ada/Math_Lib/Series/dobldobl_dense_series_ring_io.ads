with DoblDobl_Dense_Series_io;
with DoblDobl_Dense_Series_Ring;
with Abstract_Ring_io;

package DoblDobl_Dense_Series_Ring_io is
  new Abstract_Ring_io(DoblDobl_Dense_Series_Ring,
                       DoblDobl_Dense_Series_io.get,
                       DoblDobl_Dense_Series_io.put,
                       DoblDobl_Dense_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in double double precision.
