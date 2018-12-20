with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Generic_Monomials;

package Standard_Complex_Monomials is
  new Generic_Monomials(Standard_Complex_Ring,
                        Standard_Complex_Vectors,
                        Standard_Complex_VecVecs);

-- DESCRIPTION :
--   Defines monomials over the ring of standard complex numbers.
