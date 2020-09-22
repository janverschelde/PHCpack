with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_Matrices;
with DecaDobl_Complex_VecMats;
with Generic_Speelpenning_Convolutions;

package DecaDobl_Speelpenning_Convolutions is
  new Generic_Speelpenning_Convolutions(DecaDobl_Complex_Ring,
                                        DecaDobl_Complex_Vectors,
                                        DecaDobl_Complex_VecVecs,
                                        DecaDobl_Complex_Matrices,
                                        DecaDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines the evaluation and differentiation of products of power series
--   with convolutions computed over the ring of complex numbers,
--   in deca double precision.
