with HexaDobl_Complex_Numbers_io;
with HexaDobl_Complex_Ring;
with Abstract_Ring_io;

package HexaDobl_Complex_Ring_io is
  new Abstract_Ring_io(HexaDobl_Complex_Ring,
                       HexaDobl_Complex_Numbers_io.get,
                       HexaDobl_Complex_Numbers_io.put,
                       HexaDobl_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for hexa double complex numbers.
