with Multprec_Floating64_Numbers_io;
with Multprec_Floating64_Ring;
with Abstract_Ring_io;

package Multprec_Floating64_Ring_io is
  new Abstract_Ring_io(Multprec_Floating64_Ring,
                       Multprec_Floating64_Numbers_io.get,
                       Multprec_Floating64_Numbers_io.put,
                       Multprec_Floating64_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for multi-precision floating-point numbers.
