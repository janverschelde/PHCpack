with TripDobl_Complex_Ring_io;
with TripDobl_Complex_Vectors;
with Generic_Vectors_io;

package TripDobl_Complex_Vectors_io is 
  new Generic_Vectors_io(TripDobl_Complex_Ring_io,TripDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of triple double complex numbers.
