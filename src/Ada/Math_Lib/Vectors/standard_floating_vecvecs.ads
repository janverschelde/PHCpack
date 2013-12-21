with Generic_VecVecs;
with Standard_Floating_Ring;
with Standard_Floating_Vectors;

package Standard_Floating_VecVecs is 
  new Generic_VecVecs(Standard_Floating_Ring,Standard_Floating_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of standard floating numbers.
