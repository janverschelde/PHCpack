with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Generic_MatVecs;

package Double_Complex_MatVecs is 
  new Generic_MatVecs(Standard_Complex_Ring,Standard_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices of vectors over the ring of standard complex numbers.
