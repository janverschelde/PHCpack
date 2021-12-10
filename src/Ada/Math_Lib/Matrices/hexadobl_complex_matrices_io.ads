with HexaDobl_Complex_Ring_io;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_Matrices;
with Generic_Matrices_io;

package HexaDobl_Complex_Matrices_io is 
  new Generic_Matrices_io(HexaDobl_Complex_Ring_io,
                          HexaDobl_Complex_Vectors,
                          HexaDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of hexa double complex numbers.
