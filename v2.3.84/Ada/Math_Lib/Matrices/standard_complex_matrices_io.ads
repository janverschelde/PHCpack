with Standard_Complex_Ring_io;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Generic_Matrices_io;

package Standard_Complex_Matrices_io is 
  new Generic_Matrices_io(Standard_Complex_Ring_io,
                          Standard_Complex_Vectors,
                          Standard_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of standard complex numbers.
