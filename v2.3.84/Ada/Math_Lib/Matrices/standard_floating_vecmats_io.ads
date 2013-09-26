with Standard_Floating_Ring_io;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;
with Standard_Floating_VecMats;
with Generic_VecMats_io;

package Standard_Floating_VecMats_io is 
  new Generic_VecMats_io(Standard_Floating_Ring_io,
                         Standard_Floating_Vectors,
                         Standard_Floating_Matrices,
                         Standard_Floating_Matrices_io,
                         Standard_Floating_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of standard floating numbers.
