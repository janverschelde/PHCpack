with Standard_Integer64_Ring_io;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;
with Standard_Integer64_VecMats;
with Generic_VecMats_io;

package Standard_Integer64_VecMats_io is 
  new Generic_VecMats_io(Standard_Integer64_Ring_io,
                         Standard_Integer64_Vectors,
                         Standard_Integer64_Matrices,
                         Standard_Integer64_Matrices_io,
                         Standard_Integer64_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of standard integer numbers,
--   of type long long integer.
