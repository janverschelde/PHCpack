with HexaDobl_Complex_Series_Ring;
with HexaDobl_Complex_Series_Vectors;
with Generic_Matrices;

package HexaDobl_Complex_Series_Matrices is 
  new Generic_Matrices(HexaDobl_Complex_Series_Ring,
                       HexaDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of hexa double complex series.
