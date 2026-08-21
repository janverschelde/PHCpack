with Standard_Floating_Ring;
with Standard_Floating_Vectors;
with Generic_MatVecs;

package Double_Real_MatVecs is 
  new Generic_MatVecs(Standard_Floating_Ring,Standard_Floating_Vectors);

-- DESCRIPTION :
--   Defines matrices of vectors over the ring of floating-point doubles.
