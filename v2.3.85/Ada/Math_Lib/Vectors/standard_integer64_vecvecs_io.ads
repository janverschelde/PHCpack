with Standard_Integer64_Ring_io;
with Standard_Integer64_Vectors;
with Standard_Integer64_Vectors_io;
with Standard_Integer64_VecVecs;
with Generic_VecVecs_io;

package Standard_Integer64_VecVecs_io is 
  new Generic_VecVecs_io(Standard_Integer64_Ring_io,
                         Standard_Integer64_Vectors,
                         Standard_Integer64_Vectors_io,
                         Standard_Integer64_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of standard integer numbers,
--   of type long long integer.
