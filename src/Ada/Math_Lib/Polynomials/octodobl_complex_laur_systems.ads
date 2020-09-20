with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package OctoDobl_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(OctoDobl_Complex_Ring,
                                OctoDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over octo double complex numbers.
