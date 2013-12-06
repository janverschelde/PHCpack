with Standard_Integer64_Ring_io;
with Standard_Integer64_Vectors;
with Standard_Integer64_Vectors_io;
with Standard_Integer64_VecVecs;
with Lists_of_Integer64_Vectors;
with Generic_Lists_of_Vectors_io;

package Lists_of_Integer64_Vectors_io is
  new Generic_Lists_of_Vectors_io(Standard_Integer64_Ring_io,
                                  Standard_Integer64_Vectors,
                                  Standard_Integer64_Vectors_io,
                                  Standard_Integer64_VecVecs,
                                  Lists_of_Integer64_Vectors);

-- DESCRIPTION :
--   Defines input/output for lists of links to standard integer vectors,
--   with as coordinates 64-bit, or long long integers.
