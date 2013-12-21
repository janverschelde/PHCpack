with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Abstract_Ring;

package Multprec_Integer_Ring is
  new Abstract_Ring(Integer_Number,
                    Multprec_Integer_Numbers.Create(integer(0)),
                    Multprec_Integer_Numbers.Create(integer(1)),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of multi-precision integer numbers.
