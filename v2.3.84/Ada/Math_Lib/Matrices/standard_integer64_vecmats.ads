with Generic_VecMats;
with Standard_Integer64_Ring;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;

package Standard_Integer64_VecMats is 
  new Generic_VecMats(Standard_Integer64_Ring,
                      Standard_Integer64_Vectors,
                      Standard_Integer64_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of standard integer numbers,
--   of type long long integer.
