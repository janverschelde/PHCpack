with Standard_Floating_Ring_io;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Generic_Matrices_io;

package Standard_Floating_Matrices_io is 
  new Generic_Matrices_io(Standard_Floating_Ring_io,
                          Standard_Floating_Vectors,
                          Standard_Floating_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of standard floating-point numbers.
