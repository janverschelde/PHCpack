with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;

package Standard_Integer64_Ring is
  new Abstract_Ring(integer64,0,1,Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of standard integer numbers.  These are of type
--   long_long_integer, up to 2^64.
