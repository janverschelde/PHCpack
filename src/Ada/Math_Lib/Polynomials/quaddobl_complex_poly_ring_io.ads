with QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Ring;
with Abstract_Ring_io;

package QuadDobl_Complex_Poly_Ring_io is
  new Abstract_Ring_io(QuadDobl_Complex_Poly_Ring,
                       QuadDobl_Complex_Polynomials_io.get,
                       QuadDobl_Complex_Polynomials_io.put,
                       QuadDobl_Complex_Polynomials_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of complex polynomials,
--   with coefficients in quad double precision.
