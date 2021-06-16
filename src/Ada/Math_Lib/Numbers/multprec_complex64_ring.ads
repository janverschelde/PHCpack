with Multprec_Complex64_Numbers;         use Multprec_Complex64_Numbers;
with Abstract_Ring;

package Multprec_Complex64_Ring is
  new Abstract_Ring(Complex_Number,Create(integer(0)),Create(integer(1)),
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of multi-precision complex numbers.
