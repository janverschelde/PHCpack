with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_Matrices;
with OctoDobl_Complex_VecMats;
with Generic_Speelpenning_Convolutions;

package OctoDobl_Speelpenning_Convolutions is
  new Generic_Speelpenning_Convolutions(OctoDobl_Complex_Ring,
                                        OctoDobl_Complex_Vectors,
                                        OctoDobl_Complex_VecVecs,
                                        OctoDobl_Complex_Matrices,
                                        OctoDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines the evaluation and differentiation of products of power series
--   with convolutions computed over the ring of complex numbers,
--   in octo double precision.
