with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Interval_Numbers;         use Standard_Interval_Numbers;
with Abstract_Ring;

package Standard_Interval_Ring is
  new Abstract_Ring(Interval,Create(natural32(0)),Create(natural32(1)),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of standard interval numbers.
