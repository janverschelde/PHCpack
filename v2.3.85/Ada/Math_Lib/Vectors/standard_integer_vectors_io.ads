with Standard_Integer_Ring_io;
with Standard_Integer_Vectors;
with Generic_Vectors_io;

package Standard_Integer_Vectors_io is 
  new Generic_Vectors_io(Standard_Integer_Ring_io,Standard_Integer_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors over the ring of standard integer numbers.
