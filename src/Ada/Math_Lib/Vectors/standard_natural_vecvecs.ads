with Generic_VecVecs;
with Standard_Natural_Ring;
with Standard_Natural_Vectors;

package Standard_Natural_VecVecs is 
  new Generic_VecVecs(Standard_Natural_Ring,Standard_Natural_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of standard natural numbers.
