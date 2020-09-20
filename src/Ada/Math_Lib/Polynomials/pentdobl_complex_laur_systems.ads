with PentDobl_Complex_Ring;
with PentDobl_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package PentDobl_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(PentDobl_Complex_Ring,
                                PentDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over penta double complex numbers.
