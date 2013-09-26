with Standard_Integer64_Ring;
with Standard_Integer64_Vectors;
with Standard_Integer64_VecVecs;
with Generic_Lists_of_Vectors;

package Lists_of_Integer64_Vectors is 
  new Generic_Lists_of_Vectors(Standard_Integer64_Ring,
                               Standard_Integer64_Vectors,
                               Standard_Integer64_VecVecs);

-- DESCRIPTION :
--   Defines lists of links to vectors of standard integer numbers,
--   of type long long integer for 64-bit arithmetic.
