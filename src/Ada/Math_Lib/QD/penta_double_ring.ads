with Penta_Double_Numbers;                use Penta_Double_Numbers;
with Abstract_Ring;

package Penta_Double_Ring is
  new Abstract_Ring(penta_double,create(0.0),create(1.0),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of penta double numbers.
