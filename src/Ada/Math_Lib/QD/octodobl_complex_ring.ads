with OctoDobl_Complex_Numbers;           use OctoDobl_Complex_Numbers;
with Abstract_Ring;

package OctoDobl_Complex_Ring is
  new Abstract_Ring(Complex_Number,Create(integer(0)),Create(integer(1)),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of octo double complex numbers.
