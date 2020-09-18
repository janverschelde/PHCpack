with TripDobl_Complex_Series_Ring;
with TripDobl_Complex_Series_Vectors;
with Generic_Matrices;

package TripDobl_Complex_Series_Matrices is 
  new Generic_Matrices(TripDobl_Complex_Series_Ring,
                       TripDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of triple double complex series.
