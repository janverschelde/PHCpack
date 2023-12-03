with HexaDobl_Complex_Series_io;
with HexaDobl_Complex_Series_Ring;
with Abstract_Ring_io;

package HexaDobl_Complex_Series_Ring_io is
  new Abstract_Ring_io(HexaDobl_Complex_Series_Ring,
                       HexaDobl_Complex_Series_io.get,
                       HexaDobl_Complex_Series_io.put,
                       HexaDobl_Complex_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in hexa double precision.
