with DecaDobl_Complex_Numbers;           use DecaDobl_Complex_Numbers;
with Abstract_Ring;

package DecaDobl_Complex_Ring is
  new Abstract_Ring(Complex_Number,Create(integer(0)),Create(integer(1)),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of deca double complex numbers.
