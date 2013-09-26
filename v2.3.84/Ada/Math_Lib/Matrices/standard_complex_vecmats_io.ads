with Standard_Complex_Ring_io;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;
with Standard_Complex_VecMats;
with Generic_VecMats_io;

package Standard_Complex_VecMats_io is 
  new Generic_VecMats_io(Standard_Complex_Ring_io,
                         Standard_Complex_Vectors,
                         Standard_Complex_Matrices,
                         Standard_Complex_Matrices_io,
                         Standard_Complex_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of standard complex numbers.
