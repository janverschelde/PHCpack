with DecaDobl_Complex_Numbers_io;
with DecaDobl_Complex_Ring;
with Abstract_Ring_io;

package DecaDobl_Complex_Ring_io is
  new Abstract_Ring_io(DecaDobl_Complex_Ring,
                       DecaDobl_Complex_Numbers_io.get,
                       DecaDobl_Complex_Numbers_io.put,
                       DecaDobl_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for deca double complex numbers.
