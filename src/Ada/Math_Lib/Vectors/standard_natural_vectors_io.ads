with Standard_Natural_Ring_io;
with Standard_Natural_Vectors;
with Generic_Vectors_io;

package Standard_Natural_Vectors_io is 
  new Generic_Vectors_io(Standard_Natural_Ring_io,Standard_Natural_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of standard natural numbers.
