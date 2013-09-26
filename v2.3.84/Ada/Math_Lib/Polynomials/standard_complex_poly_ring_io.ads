with Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Ring;
with Abstract_Ring_io;

package Standard_Complex_Poly_Ring_io is
  new Abstract_Ring_io(Standard_Complex_Poly_Ring,
                       Standard_Complex_Polynomials_io.get,
                       Standard_Complex_Polynomials_io.put,
                       Standard_Complex_Polynomials_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of standard complex polynomials.
