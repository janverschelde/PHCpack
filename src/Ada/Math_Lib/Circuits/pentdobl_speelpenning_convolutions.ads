with PentDobl_Complex_Ring;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_Matrices;
with PentDobl_Complex_VecMats;
with Generic_Speelpenning_Convolutions;

package PentDobl_Speelpenning_Convolutions is
  new Generic_Speelpenning_Convolutions(PentDobl_Complex_Ring,
                                        PentDobl_Complex_Vectors,
                                        PentDobl_Complex_VecVecs,
                                        PentDobl_Complex_Matrices,
                                        PentDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines the evaluation and differentiation of products of power series
--   with convolutions computed over the ring of complex numbers,
--   in penta double precision.
