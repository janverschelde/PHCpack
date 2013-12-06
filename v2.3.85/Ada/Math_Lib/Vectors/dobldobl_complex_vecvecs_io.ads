with DoblDobl_Complex_Ring_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with Generic_VecVecs_io;

package DoblDobl_Complex_VecVecs_io is 
  new Generic_VecVecs_io(DoblDobl_Complex_Ring_io,
                         DoblDobl_Complex_Vectors,
                         DoblDobl_Complex_Vectors_io,
                         DoblDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors
--   of double double complex numbers.
