with DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Ring;
with Abstract_Ring_io;

package DoblDobl_Complex_Ring_io is
  new Abstract_Ring_io(DoblDobl_Complex_Ring,
                       DoblDobl_Complex_Numbers_io.get,
                       DoblDobl_Complex_Numbers_io.put,
                       DoblDobl_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for double double complex numbers.
