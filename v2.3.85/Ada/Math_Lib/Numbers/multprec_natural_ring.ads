with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Abstract_Ring;

package Multprec_Natural_Ring is
  new Abstract_Ring(Natural_Number,
                    Multprec_Natural_Numbers.Create(integer(0)),
                    Multprec_Natural_Numbers.Create(integer(1)),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of multi-precision natural numbers.
