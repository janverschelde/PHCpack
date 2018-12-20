with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Monomials;
with Generic_Monomial_Vectors;

package DoblDobl_Monomial_Vectors is
  new Generic_Monomial_Vectors(DoblDobl_Complex_Ring,
                               DoblDobl_Complex_Vectors,
                               DoblDobl_Complex_VecVecs,
                               DoblDobl_Complex_Monomials);

-- DESCRIPTION :
--   Defines monomial vectors over the ring of double double complex numbers.
