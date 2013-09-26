with Standard_Complex_Ring_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Lists_of_Complex_Vectors;
with Generic_Lists_of_Vectors_io;

package Lists_of_Complex_Vectors_io is
  new Generic_Lists_of_Vectors_io(Standard_Complex_Ring_io,
                                  Standard_Complex_Vectors,
                                  Standard_Complex_Vectors_io,
                                  Standard_Complex_VecVecs,
                                  Lists_of_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output for lists of links to standard complex vectors.
