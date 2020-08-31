with TripDobl_Complex_Ring_io;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Vectors_io;
with TripDobl_Complex_VecVecs;
with Generic_VecVecs_io;

package TripDobl_Complex_VecVecs_io is 
  new Generic_VecVecs_io(TripDobl_Complex_Ring_io,
                         TripDobl_Complex_Vectors,
                         TripDobl_Complex_Vectors_io,
                         TripDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors
--   of triple double complex numbers.
