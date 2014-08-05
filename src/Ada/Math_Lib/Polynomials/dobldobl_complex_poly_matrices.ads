with DoblDobl_Complex_Poly_Ring;         use DoblDobl_Complex_Poly_Ring;
with DoblDobl_Complex_Poly_Vectors;
with Generic_Matrices;

package DoblDobl_Complex_Poly_Matrices is
  new Generic_Matrices(DoblDobl_Complex_Poly_Ring,
                       DoblDobl_Complex_Poly_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of double double complex polynomials.
