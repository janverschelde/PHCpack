with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_VecVecs;
with HexaDobl_Complex_Matrices;
with HexaDobl_Complex_VecMats;
with Generic_Speelpenning_Convolutions;

package HexaDobl_Speelpenning_Convolutions is
  new Generic_Speelpenning_Convolutions(HexaDobl_Complex_Ring,
                                        HexaDobl_Complex_Vectors,
                                        HexaDobl_Complex_VecVecs,
                                        HexaDobl_Complex_Matrices,
                                        HexaDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines the evaluation and differentiation of products of power series
--   with convolutions computed over the ring of complex numbers,
--   in hexa double precision.
