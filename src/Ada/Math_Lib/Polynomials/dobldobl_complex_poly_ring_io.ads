with DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Ring;
with Abstract_Ring_io;

package DoblDobl_Complex_Poly_Ring_io is
  new Abstract_Ring_io(DoblDobl_Complex_Poly_Ring,
                       DoblDobl_Complex_Polynomials_io.get,
                       DoblDobl_Complex_Polynomials_io.put,
                       DoblDobl_Complex_Polynomials_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of complex polynomials,
--   with coefficients in double double precision.
