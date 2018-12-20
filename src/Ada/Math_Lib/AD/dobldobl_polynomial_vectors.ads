with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Monomials;
with DoblDobl_Monomial_Vectors;
with Generic_Polynomial_Vectors;

package DoblDobl_Polynomial_Vectors is
  new Generic_Polynomial_Vectors(DoblDobl_Complex_Ring,
                                 DoblDobl_Complex_Vectors,
                                 DoblDobl_Complex_VecVecs,
                                 DoblDobl_Complex_Matrices,
                                 DoblDobl_Complex_Monomials,
                                 DoblDobl_Monomial_Vectors);

-- DESCRIPTION :
--   Defines polynomial vectors over the ring of double double complex numbers.
