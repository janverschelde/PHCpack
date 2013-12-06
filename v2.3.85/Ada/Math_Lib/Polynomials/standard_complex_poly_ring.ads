with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Abstract_Ring;

package Standard_Complex_Poly_Ring is
  new Abstract_Ring(Poly,Null_Poly,One_Poly,
                    Create,Equal,Copy,"+","+","-","-","*",
                    Add,Sub,Min,Mul,Clear);

-- DESCRIPTION :
--   Defines the ring of polynomials in several variables with as
--   coefficients the standard complex numbers.
