with Standard_Floating_Ring_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Lists_of_Floating_Vectors;
with Lists_of_Floating_Vectors_io;
with Arrays_of_Floating_Vector_Lists;
with Generic_Arrays_of_Vector_Lists_io;

package Arrays_of_Floating_Vector_Lists_io is
  new Generic_Arrays_of_Vector_Lists_io(Standard_Floating_Ring_io,
                                        Standard_Floating_Vectors,
                                        Standard_Floating_Vectors_io,
                                        Standard_Floating_VecVecs,
                                        Lists_of_Floating_Vectors,
                                        Lists_of_Floating_Vectors_io,
                                        Arrays_of_Floating_Vector_Lists);

-- DESCRIPTION :
--   Defines input/output for arrays of lists of links to floating vectors.
