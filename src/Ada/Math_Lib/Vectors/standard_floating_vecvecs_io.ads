with Standard_Floating_Ring_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Generic_VecVecs_io;

package Standard_Floating_VecVecs_io is 
  new Generic_VecVecs_io(Standard_Floating_Ring_io,
                         Standard_Floating_Vectors,
                         Standard_Floating_Vectors_io,
                         Standard_Floating_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of standard floating numbers.
