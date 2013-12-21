with Standard_Integer64_Ring_io;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;
with Generic_Matrices_io;

package Standard_Integer64_Matrices_io is 
  new Generic_Matrices_io(Standard_Integer64_Ring_io,
                          Standard_Integer64_Vectors,
                          Standard_Integer64_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of standard integer numbers,
--   of type long long integer.
