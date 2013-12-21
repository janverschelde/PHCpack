with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Abstract_Ring;

package Standard_Natural64_Ring is
  new Abstract_Ring(natural64,0,1,Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of standard natural numbers, of range 0..2^64.
