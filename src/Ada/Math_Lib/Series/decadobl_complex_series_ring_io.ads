with DecaDobl_Complex_Series_io;
with DecaDobl_Complex_Series_Ring;
with Abstract_Ring_io;

package DecaDobl_Complex_Series_Ring_io is
  new Abstract_Ring_io(DecaDobl_Complex_Series_Ring,
                       DecaDobl_Complex_Series_io.get,
                       DecaDobl_Complex_Series_io.put,
                       DecaDobl_Complex_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in deca double precision.
