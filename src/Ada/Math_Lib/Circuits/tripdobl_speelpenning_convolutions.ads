with TripDobl_Complex_Ring;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_VecMats;
with Generic_Speelpenning_Convolutions;

package TripDobl_Speelpenning_Convolutions is
  new Generic_Speelpenning_Convolutions(TripDobl_Complex_Ring,
                                        TripDobl_Complex_Vectors,
                                        TripDobl_Complex_VecVecs,
                                        TripDobl_Complex_Matrices,
                                        TripDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines the evaluation and differentiation of products of power series
--   with convolutions computed over the ring of complex numbers,
--   in triple double precision.
