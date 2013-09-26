with QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Ring;
with Abstract_Ring_io;

package QuadDobl_Complex_Ring_io is
  new Abstract_Ring_io(QuadDobl_Complex_Ring,
                       QuadDobl_Complex_Numbers_io.get,
                       QuadDobl_Complex_Numbers_io.put,
                       QuadDobl_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for quad double complex numbers.
