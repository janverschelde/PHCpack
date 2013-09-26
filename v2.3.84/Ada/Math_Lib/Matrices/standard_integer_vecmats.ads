with Generic_VecMats;
with Standard_Integer_Ring;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;

package Standard_Integer_VecMats is 
  new Generic_VecMats(Standard_Integer_Ring,
                      Standard_Integer_Vectors,
                      Standard_Integer_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of standard integer numbers.
