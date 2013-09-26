with Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Ring;
with Abstract_Ring_io;

package Multprec_Complex_Poly_Ring_io is
  new Abstract_Ring_io(Multprec_Complex_Poly_Ring,
                       Multprec_Complex_Polynomials_io.get,
                       Multprec_Complex_Polynomials_io.put,
                       Multprec_Complex_Polynomials_io.put);

-- DESCRIPTION :
--   Defines input/output for the ring of multiprecision complex polynomials.
