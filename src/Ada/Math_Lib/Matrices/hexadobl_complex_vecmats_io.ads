with HexaDobl_Complex_Ring_io;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_Matrices;
with HexaDobl_Complex_Matrices_io;
with HexaDobl_Complex_VecMats;
with Generic_VecMats_io;

package HexaDobl_Complex_VecMats_io is 
  new Generic_VecMats_io(HexaDobl_Complex_Ring_io,
                         HexaDobl_Complex_Vectors,
                         HexaDobl_Complex_Matrices,
                         HexaDobl_Complex_Matrices_io,
                         HexaDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of complex numbers
--   in hexa double precision.
