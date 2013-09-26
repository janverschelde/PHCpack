with Standard_Integer_Ring;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Generic_Lists_of_Vectors;

package Lists_of_Integer_Vectors is 
  new Generic_Lists_of_Vectors(Standard_Integer_Ring,
                               Standard_Integer_Vectors,
                               Standard_Integer_VecVecs);

-- DESCRIPTION :
--   Defines lists of links to vectors of standard integer numbers.
