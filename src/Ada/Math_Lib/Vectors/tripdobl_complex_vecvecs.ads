with Generic_VecVecs;
with TripDobl_Complex_Ring;
with TripDobl_Complex_Vectors;

package TripDobl_Complex_VecVecs is 
  new Generic_VecVecs(TripDobl_Complex_Ring,TripDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of triple double complex numbers.
