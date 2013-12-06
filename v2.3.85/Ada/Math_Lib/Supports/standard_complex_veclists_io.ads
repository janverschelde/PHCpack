with Standard_Complex_Ring_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecLists;
with Generic_Lists_of_Vectors_io;

package Standard_Complex_VecLists_io is
  new Generic_Lists_of_Vectors_io(Standard_Complex_Ring_io,
                                  Standard_Complex_Vectors,
                                  Standard_Complex_Vectors_io,
                                  Standard_Complex_VecVecs,
                                  Standard_Complex_VecLists);

-- DESCRIPTION :
--   Defines input/output for lists of standard complex vectors.
