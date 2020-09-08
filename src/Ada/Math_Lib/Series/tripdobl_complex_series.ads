with TripDobl_Complex_Numbers;            use TripDobl_Complex_Numbers;
with TripDobl_Complex_Ring;
with TripDobl_Complex_Vectors;
with Generic_Dense_Series;

package TripDobl_Complex_Series is
  new Generic_Dense_Series(TripDobl_Complex_Ring,
                           TripDobl_Complex_Vectors,"/",Div);

-- DESCRIPTION :
--   Defines series using vectors over the ring of complex numbers,
--   with the division, in triple double precision.
