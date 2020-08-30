with TripDobl_Complex_Ring_io;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Matrices;
with Generic_Matrices_io;

package TripDobl_Complex_Matrices_io is 
  new Generic_Matrices_io(TripDobl_Complex_Ring_io,
                          TripDobl_Complex_Vectors,
                          TripDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of triple double complex numbers.
