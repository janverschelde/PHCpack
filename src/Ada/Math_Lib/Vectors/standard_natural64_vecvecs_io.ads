with Standard_Natural64_Ring_io;
with Standard_Natural64_Vectors;
with Standard_Natural64_Vectors_io;
with Standard_Natural64_VecVecs;
with Generic_VecVecs_io;

package Standard_Natural64_VecVecs_io is 
  new Generic_VecVecs_io(Standard_Natural64_Ring_io,
                         Standard_Natural64_Vectors,
                         Standard_Natural64_Vectors_io,
                         Standard_Natural64_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of standard natural numbers.
