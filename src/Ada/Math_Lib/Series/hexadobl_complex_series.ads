with HexaDobl_Complex_Numbers;            use HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Vectors;
with Generic_Dense_Series;

package HexaDobl_Complex_Series is
  new Generic_Dense_Series(HexaDobl_Complex_Ring,
                           HexaDobl_Complex_Vectors,"/",Div);

-- DESCRIPTION :
--   Defines series using vectors over the ring of complex numbers,
--   with the division, in hexa double precision.
