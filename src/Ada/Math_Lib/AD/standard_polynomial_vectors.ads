with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Monomials;
with Standard_Monomial_Vectors;
with Generic_Polynomial_Vectors;

package Standard_Polynomial_Vectors is
  new Generic_Polynomial_Vectors(Standard_Complex_Ring,
                                 Standard_Complex_Vectors,
                                 Standard_Complex_VecVecs,
                                 Standard_Complex_Matrices,
                                 Standard_Complex_Monomials,
                                 Standard_Monomial_Vectors);

-- DESCRIPTION :
--   Defines polynomial vectors over the ring of standard complex numbers.
