with Standard_Natural64_Ring_io;
with Standard_Natural64_Vectors;
with Generic_Vectors_io;

package Standard_Natural64_Vectors_io is 
  new Generic_Vectors_io(Standard_Natural64_Ring_io,
                         Standard_Natural64_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of standard natural numbers.
