with QuadDobl_Complex_Poly_Ring;         use QuadDobl_Complex_Poly_Ring;
with QuadDobl_Complex_Poly_Vectors;
with Generic_Matrices;

package QuadDobl_Complex_Poly_Matrices is
  new Generic_Matrices(QuadDobl_Complex_Poly_Ring,
                       QuadDobl_Complex_Poly_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of quad double complex polynomials.
