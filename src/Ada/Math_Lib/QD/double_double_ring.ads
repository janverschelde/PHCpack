with Double_Double_Numbers;               use Double_Double_Numbers;
with Abstract_Ring;

package Double_Double_Ring is
  new Abstract_Ring(double_double,create(0.0),create(1.0),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of double double numbers.
