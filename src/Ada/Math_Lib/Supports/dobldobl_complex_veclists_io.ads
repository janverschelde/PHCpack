with DoblDobl_Complex_Ring_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecLists;
with Generic_Lists_of_Vectors_io;

package DoblDobl_Complex_VecLists_io is
  new Generic_Lists_of_Vectors_io(DoblDobl_Complex_Ring_io,
                                  DoblDobl_Complex_Vectors,
                                  DoblDobl_Complex_Vectors_io,
                                  DoblDobl_Complex_VecVecs,
                                  DoblDobl_Complex_VecLists);

-- DESCRIPTION :
--   Defines input/output for lists of double double complex vectors.
