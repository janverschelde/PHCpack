with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package DoblDobl_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(DoblDobl_Complex_Ring,
                                DoblDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over double double complex numbers.
