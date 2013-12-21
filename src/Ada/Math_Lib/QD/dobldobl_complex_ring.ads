with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Abstract_Ring;

package DoblDobl_Complex_Ring is
  new Abstract_Ring(Complex_Number,Create(integer(0)),Create(integer(1)),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of double double complex numbers.
