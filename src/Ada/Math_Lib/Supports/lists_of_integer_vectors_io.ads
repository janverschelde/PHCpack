with Standard_Integer_Ring_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;
with Generic_Lists_of_Vectors_io;

package Lists_of_Integer_Vectors_io is
  new Generic_Lists_of_Vectors_io(Standard_Integer_Ring_io,
                                  Standard_Integer_Vectors,
                                  Standard_Integer_Vectors_io,
                                  Standard_Integer_VecVecs,
                                  Lists_of_Integer_Vectors);

-- DESCRIPTION :
--   Defines input/output for lists of links to standard integer vectors.
