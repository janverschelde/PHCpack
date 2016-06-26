with QuadDobl_Dense_Series_io;
with QuadDobl_Dense_Series_Ring;
with Abstract_Ring_io;

package QuadDobl_Dense_Series_Ring_io is
  new Abstract_Ring_io(QuadDobl_Dense_Series_Ring,
                       QuadDobl_Dense_Series_io.get,
                       QuadDobl_Dense_Series_io.put,
                       QuadDobl_Dense_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in quad double precision.
