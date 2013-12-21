with Standard_Integer_Ring_io;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Generic_Matrices_io;

package Standard_Integer_Matrices_io is 
  new Generic_Matrices_io(Standard_Integer_Ring_io,
                          Standard_Integer_Vectors,
                          Standard_Integer_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of standard integer numbers.
