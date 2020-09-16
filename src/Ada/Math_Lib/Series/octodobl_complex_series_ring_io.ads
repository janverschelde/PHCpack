with OctoDobl_Complex_Series_io;
with OctoDobl_Complex_Series_Ring;
with Abstract_Ring_io;

package OctoDobl_Complex_Series_Ring_io is
  new Abstract_Ring_io(OctoDobl_Complex_Series_Ring,
                       OctoDobl_Complex_Series_io.get,
                       OctoDobl_Complex_Series_io.put,
                       OctoDobl_Complex_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in octo double precision.
