with Generic_VecVecs;
with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Vectors;

package HexaDobl_Complex_VecVecs is 
  new Generic_VecVecs(HexaDobl_Complex_Ring,HexaDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of hexa double complex numbers.
