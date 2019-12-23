with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with Generic_Speelpenning_Convolutions;

package QuadDobl_Speelpenning_Convolutions is
  new Generic_Speelpenning_Convolutions(QuadDobl_Complex_Ring,
                                        QuadDobl_Complex_Vectors,
                                        QuadDobl_Complex_VecVecs,
                                        QuadDobl_Complex_Matrices,
                                        QuadDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines the evaluation and differentiation of products of power series
--   with convolutions computed over the ring of complex numbers,
--   in quad double precision.
