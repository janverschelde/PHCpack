with TripDobl_Complex_Ring_io;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_Matrices_io;
with TripDobl_Complex_VecMats;
with Generic_VecMats_io;

package TripDobl_Complex_VecMats_io is 
  new Generic_VecMats_io(TripDobl_Complex_Ring_io,
                         TripDobl_Complex_Vectors,
                         TripDobl_Complex_Matrices,
                         TripDobl_Complex_Matrices_io,
                         TripDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of complex numbers
--   in triple double precision.
