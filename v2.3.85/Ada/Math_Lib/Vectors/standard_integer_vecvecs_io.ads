with Standard_Integer_Ring_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Generic_VecVecs_io;

package Standard_Integer_VecVecs_io is 
  new Generic_VecVecs_io(Standard_Integer_Ring_io,
                         Standard_Integer_Vectors,
                         Standard_Integer_Vectors_io,
                         Standard_Integer_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of standard integer numbers.
