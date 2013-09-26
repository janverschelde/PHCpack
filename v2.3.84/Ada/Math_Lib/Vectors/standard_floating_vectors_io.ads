with Standard_Floating_Ring_io;
with Standard_Floating_Vectors;
with Generic_Vectors_io;

package Standard_Floating_Vectors_io is 
  new Generic_Vectors_io(Standard_Floating_Ring_io,Standard_Floating_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of standard floating-point numbers.
