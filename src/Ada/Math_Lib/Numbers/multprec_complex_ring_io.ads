with Multprec_Complex_Numbers_io;
with Multprec_Complex_Ring;
with Abstract_Ring_io;

package Multprec_Complex_Ring_io is
  new Abstract_Ring_io(Multprec_Complex_Ring,
                       Multprec_Complex_Numbers_io.get,
                       Multprec_Complex_Numbers_io.put,
                       Multprec_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for multi-precision complex numbers.
