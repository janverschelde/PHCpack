with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package HexaDobl_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(HexaDobl_Complex_Ring,
                                HexaDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over hexa double complex numbers.
