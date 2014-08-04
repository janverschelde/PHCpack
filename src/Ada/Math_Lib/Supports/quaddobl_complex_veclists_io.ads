with QuadDobl_Complex_Ring_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecLists;
with Generic_Lists_of_Vectors_io;

package QuadDobl_Complex_VecLists_io is
  new Generic_Lists_of_Vectors_io(QuadDobl_Complex_Ring_io,
                                  QuadDobl_Complex_Vectors,
                                  QuadDobl_Complex_Vectors_io,
                                  QuadDobl_Complex_VecVecs,
                                  QuadDobl_Complex_VecLists);

-- DESCRIPTION :
--   Defines input/output for lists of quad double complex vectors.
