with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Monomials;
with Generic_Monomial_Vectors;

package Standard_Monomial_Vectors is
  new Generic_Monomial_Vectors(Standard_Complex_Ring,
                               Standard_Complex_Vectors,
                               Standard_Complex_VecVecs,
                               Standard_Complex_Monomials);

-- DESCRIPTION :
--   Defines monomial vectors over the ring of standard complex numbers.
