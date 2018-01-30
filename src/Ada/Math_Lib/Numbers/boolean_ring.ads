with Boolean_Numbers;                    use Boolean_Numbers;
with Abstract_Ring;

package Boolean_Ring is
  new Abstract_Ring(boolean,false,true,Create,Equal,Copy,
                    "+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of Boolean numbers.
