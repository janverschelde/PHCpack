with Standard_Natural64_Ring;
with Generic_Vectors;

package Standard_Natural64_Vectors is 
  new Generic_Vectors(Standard_Natural64_Ring);

-- DESCRIPTION :
--   Defines vectors over the ring of standard natural numbers,
--   using 64-bit integer arithmetic.
