with Standard_Complex_Ring_io;
with Standard_Complex_Vectors;
with Generic_Vectors_io;

package Standard_Complex_Vectors_io is 
  new Generic_Vectors_io(Standard_Complex_Ring_io,Standard_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of standard complex numbers.
