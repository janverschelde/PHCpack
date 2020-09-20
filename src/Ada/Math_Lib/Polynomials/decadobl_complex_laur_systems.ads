with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package DecaDobl_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(DecaDobl_Complex_Ring,
                                DecaDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over deca double complex numbers.
