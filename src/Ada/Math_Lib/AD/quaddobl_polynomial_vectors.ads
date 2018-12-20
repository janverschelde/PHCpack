with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Monomials;
with QuadDobl_Monomial_Vectors;
with Generic_Polynomial_Vectors;

package QuadDobl_Polynomial_Vectors is
  new Generic_Polynomial_Vectors(QuadDobl_Complex_Ring,
                                 QuadDobl_Complex_Vectors,
                                 QuadDobl_Complex_VecVecs,
                                 QuadDobl_Complex_Matrices,
                                 QuadDobl_Complex_Monomials,
                                 QuadDobl_Monomial_Vectors);

-- DESCRIPTION :
--   Defines polynomial vectors over the ring of quad double complex numbers.
