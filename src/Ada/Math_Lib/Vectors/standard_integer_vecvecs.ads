with Generic_VecVecs;
with Standard_Integer_Ring;
with Standard_Integer_Vectors;

package Standard_Integer_VecVecs is 
  new Generic_VecVecs(Standard_Integer_Ring,Standard_Integer_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of standard integer numbers.
