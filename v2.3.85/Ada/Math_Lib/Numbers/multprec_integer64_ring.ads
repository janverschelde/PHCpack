with Multprec_Integer64_Numbers;         use Multprec_Integer64_Numbers;
with Abstract_Ring;

package Multprec_Integer64_Ring is
  new Abstract_Ring(Integer_Number,
                    Multprec_Integer64_Numbers.Create64(0),
                    Multprec_Integer64_Numbers.Create64(1),
                    Multprec_Integer64_Numbers.Create,
                    Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of multi-precision integer numbers,
--   using standard 64-bit integer arithmetic.
