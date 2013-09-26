with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Abstract_Ring;

package Standard_Floating_Ring is
  new Abstract_Ring(double_float,0.0,1.0,
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of standard double floating-point numbers.
