with Standard_Integer_Ring;              use Standard_Integer_Ring;
with Standard_Integer_Vectors;
with Generic_Matrices;

package Standard_Integer_Matrices is
  new Generic_Matrices(Standard_Integer_Ring,
                       Standard_Integer_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard integer numbers.
