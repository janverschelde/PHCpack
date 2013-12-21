with Standard_Integer64_Ring_io;
with Standard_Integer64_Vectors;
with Generic_Vectors_io;

package Standard_Integer64_Vectors_io is 
  new Generic_Vectors_io(Standard_Integer64_Ring_io,
                         Standard_Integer64_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors over the ring of standard integer numbers,
--   of type long long integer.
