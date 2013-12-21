with Multprec_Floating_Numbers_io;
with Multprec_Floating_Ring;
with Abstract_Ring_io;

package Multprec_Floating_Ring_io is
  new Abstract_Ring_io(Multprec_Floating_Ring,
                       Multprec_Floating_Numbers_io.get,
                       Multprec_Floating_Numbers_io.put,
                       Multprec_Floating_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for multi-precision floating-point numbers.
