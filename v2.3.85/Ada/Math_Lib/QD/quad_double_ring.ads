with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Abstract_Ring;

package Quad_Double_Ring is
  new Abstract_Ring(quad_double,create(0.0),create(1.0),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of quad double numbers.
