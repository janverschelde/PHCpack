with TripDobl_Complex_Ring;              use TripDobl_Complex_Ring;
with TripDobl_Complex_Vectors;
with Generic_Matrices;

package TripDobl_Complex_Matrices is
  new Generic_Matrices(TripDobl_Complex_Ring,TripDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of triple double complex numbers.
