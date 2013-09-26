with Standard_Integer_Ring_io;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;
with Standard_Integer_VecMats;
with Generic_VecMats_io;

package Standard_Integer_VecMats_io is 
  new Generic_VecMats_io(Standard_Integer_Ring_io,
                         Standard_Integer_Vectors,
                         Standard_Integer_Matrices,
                         Standard_Integer_Matrices_io,
                         Standard_Integer_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of standard integer numbers.
