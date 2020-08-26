with Octo_Double_Numbers;                 use Octo_Double_Numbers;
with Abstract_Ring;

package Octo_Double_Ring is
  new Abstract_Ring(octo_double,create(0.0),create(1.0),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of octo double numbers.
