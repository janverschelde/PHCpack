with Hexa_Double_Numbers;                 use Hexa_Double_Numbers;
with Abstract_Ring;

package Hexa_Double_Ring is
  new Abstract_Ring(hexa_double,create(0.0),create(1.0),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of hexa double numbers.
