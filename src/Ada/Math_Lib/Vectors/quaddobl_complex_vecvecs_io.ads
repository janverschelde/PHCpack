with QuadDobl_Complex_Ring_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with Generic_VecVecs_io;

package QuadDobl_Complex_VecVecs_io is 
  new Generic_VecVecs_io(QuadDobl_Complex_Ring_io,
                         QuadDobl_Complex_Vectors,
                         QuadDobl_Complex_Vectors_io,
                         QuadDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors
--   of quad double complex numbers.
