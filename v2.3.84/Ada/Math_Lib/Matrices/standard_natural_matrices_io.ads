with Standard_Natural_Ring_io;
with Standard_Natural_Vectors;
with Standard_Natural_Matrices;
with Generic_Matrices_io;

package Standard_Natural_Matrices_io is 
  new Generic_Matrices_io(Standard_Natural_Ring_io,
                          Standard_Natural_Vectors,
                          Standard_Natural_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of standard natural numbers.
