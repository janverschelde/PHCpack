with Standard_Natural_Ring_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Generic_VecVecs_io;

package Standard_Natural_VecVecs_io is 
  new Generic_VecVecs_io(Standard_Natural_Ring_io,
                         Standard_Natural_Vectors,
                         Standard_Natural_Vectors_io,
                         Standard_Natural_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of standard natural numbers.
