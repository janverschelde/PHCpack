with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Monomials;
with Generic_Monomial_Vectors;

package QuadDobl_Monomial_Vectors is
  new Generic_Monomial_Vectors(QuadDobl_Complex_Ring,
                               QuadDobl_Complex_Vectors,
                               QuadDobl_Complex_VecVecs,
                               QuadDobl_Complex_Monomials);

-- DESCRIPTION :
--   Defines monomial vectors over the ring of quad double complex numbers.
