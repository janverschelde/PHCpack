with DecaDobl_Complex_Ring_io;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Vectors_io;
with DecaDobl_Complex_VecVecs;
with Generic_VecVecs_io;

package DecaDobl_Complex_VecVecs_io is 
  new Generic_VecVecs_io(DecaDobl_Complex_Ring_io,
                         DecaDobl_Complex_Vectors,
                         DecaDobl_Complex_Vectors_io,
                         DecaDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors
--   of deca double complex numbers.
