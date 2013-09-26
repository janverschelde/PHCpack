with Multprec_Floating_Numbers;           use Multprec_Floating_Numbers;
with Abstract_Ring;

package Multprec_Floating_Ring is
  new Abstract_Ring(Floating_Number,Create(0.0),Create(1.0),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of multi-precision floating-point numbers.
