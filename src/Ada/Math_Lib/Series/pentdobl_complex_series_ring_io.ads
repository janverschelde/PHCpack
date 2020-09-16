with PentDobl_Complex_Series_io;
with PentDobl_Complex_Series_Ring;
with Abstract_Ring_io;

package PentDobl_Complex_Series_Ring_io is
  new Abstract_Ring_io(PentDobl_Complex_Series_Ring,
                       PentDobl_Complex_Series_io.get,
                       PentDobl_Complex_Series_io.put,
                       PentDobl_Complex_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in penta double precision.
