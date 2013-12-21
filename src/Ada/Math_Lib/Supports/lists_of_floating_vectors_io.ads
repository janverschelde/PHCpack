with Standard_Floating_Ring_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Lists_of_Floating_Vectors;
with Generic_Lists_of_Vectors_io;

package Lists_of_Floating_Vectors_io is
  new Generic_Lists_of_Vectors_io(Standard_Floating_Ring_io,
                                  Standard_Floating_Vectors,
                                  Standard_Floating_Vectors_io,
                                  Standard_Floating_VecVecs,
                                  Lists_of_Floating_Vectors);

-- DESCRIPTION :
--   Defines input/output for lists of links to standard floating-point vectors.
