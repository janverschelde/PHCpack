with HexaDobl_Complex_Ring_io;
with HexaDobl_Complex_Vectors;
with Generic_Vectors_io;

package HexaDobl_Complex_Vectors_io is 
  new Generic_Vectors_io(HexaDobl_Complex_Ring_io,HexaDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of hexa double complex numbers.
