with OctoDobl_Complex_Ring_io;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Vectors_io;
with OctoDobl_Complex_VecVecs;
with Generic_VecVecs_io;

package OctoDobl_Complex_VecVecs_io is 
  new Generic_VecVecs_io(OctoDobl_Complex_Ring_io,
                         OctoDobl_Complex_Vectors,
                         OctoDobl_Complex_Vectors_io,
                         OctoDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors
--   of octo double complex numbers.
