with HexaDobl_Complex_Ring;              use HexaDobl_Complex_Ring;
with HexaDobl_Complex_Vectors;
with Generic_Matrices;

package HexaDobl_Complex_Matrices is
  new Generic_Matrices(HexaDobl_Complex_Ring,HexaDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of hexa double complex numbers.
