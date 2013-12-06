with Multprec_Natural_Numbers_io;
with Multprec_Natural_Ring;
with Abstract_Ring_io;

package Multprec_Natural_Ring_io is
  new Abstract_Ring_io(Multprec_Natural_Ring,
                       Multprec_Natural_Numbers_io.get,
                       Multprec_Natural_Numbers_io.put,
                       Multprec_Natural_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of multi-precision natural numbers.
