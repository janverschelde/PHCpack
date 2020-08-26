with Triple_Double_Numbers;               use Triple_Double_Numbers;
with Abstract_Ring;

package Triple_Double_Ring is
  new Abstract_Ring(triple_double,create(0.0),create(1.0),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of triple double numbers.
