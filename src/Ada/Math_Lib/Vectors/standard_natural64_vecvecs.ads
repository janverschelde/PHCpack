with Generic_VecVecs;
with Standard_Natural64_Ring;
with Standard_Natural64_Vectors;

package Standard_Natural64_VecVecs is 
  new Generic_VecVecs(Standard_Natural64_Ring,Standard_Natural64_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of standard natural numbers,
--   using 64-bit integer arithmetic.
