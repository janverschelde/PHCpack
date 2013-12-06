with Multprec_Integer_Numbers_io;
with Multprec_Integer_Ring;
with Abstract_Ring_io;

package Multprec_Integer_Ring_io is
  new Abstract_Ring_io(Multprec_Integer_Ring,
                       Multprec_Integer_Numbers_io.get,
                       Multprec_Integer_Numbers_io.put,
                       Multprec_Integer_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of multi-precision integer numbers.
