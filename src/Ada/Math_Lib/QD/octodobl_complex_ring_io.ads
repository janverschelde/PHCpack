with OctoDobl_Complex_Numbers_io;
with OctoDobl_Complex_Ring;
with Abstract_Ring_io;

package OctoDobl_Complex_Ring_io is
  new Abstract_Ring_io(OctoDobl_Complex_Ring,
                       OctoDobl_Complex_Numbers_io.get,
                       OctoDobl_Complex_Numbers_io.put,
                       OctoDobl_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for octo double complex numbers.
