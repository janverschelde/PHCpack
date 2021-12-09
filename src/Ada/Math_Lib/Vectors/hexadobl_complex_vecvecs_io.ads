with HexaDobl_Complex_Ring_io;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_Vectors_io;
with HexaDobl_Complex_VecVecs;
with Generic_VecVecs_io;

package HexaDobl_Complex_VecVecs_io is 
  new Generic_VecVecs_io(HexaDobl_Complex_Ring_io,
                         HexaDobl_Complex_Vectors,
                         HexaDobl_Complex_Vectors_io,
                         HexaDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors
--   of hexa double complex numbers.
