with PentDobl_Complex_Numbers_io;
with PentDobl_Complex_Ring;
with Abstract_Ring_io;

package PentDobl_Complex_Ring_io is
  new Abstract_Ring_io(PentDobl_Complex_Ring,
                       PentDobl_Complex_Numbers_io.get,
                       PentDobl_Complex_Numbers_io.put,
                       PentDobl_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for penta double complex numbers.
