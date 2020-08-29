with Deca_Double_Numbers;                 use Deca_Double_Numbers;
with Abstract_Ring;

package Deca_Double_Ring is
  new Abstract_Ring(deca_double,create(0.0),create(1.0),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of deca double numbers.
