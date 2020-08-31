with PentDobl_Complex_Ring_io;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Vectors_io;
with PentDobl_Complex_VecVecs;
with Generic_VecVecs_io;

package PentDobl_Complex_VecVecs_io is 
  new Generic_VecVecs_io(PentDobl_Complex_Ring_io,
                         PentDobl_Complex_Vectors,
                         PentDobl_Complex_Vectors_io,
                         PentDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors
--   of penta double complex numbers.
