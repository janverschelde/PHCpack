with Generic_VecVecs;
with Standard_Integer64_Ring;
with Standard_Integer64_Vectors;

package Standard_Integer64_VecVecs is 
  new Generic_VecVecs(Standard_Integer64_Ring,Standard_Integer64_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of standard integer numbers,
--   of type long long integer.
