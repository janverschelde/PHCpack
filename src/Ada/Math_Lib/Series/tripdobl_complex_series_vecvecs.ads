with Generic_VecVecs;
with TripDobl_Complex_Series_Ring;
with TripDobl_Complex_Series_Vectors;

package TripDobl_Complex_Series_VecVecs is 
  new Generic_VecVecs(TripDobl_Complex_Series_Ring,
                      TripDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with triple double complex numbers as coefficients.
