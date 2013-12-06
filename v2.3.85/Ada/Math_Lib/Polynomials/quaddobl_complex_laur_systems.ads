with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package QuadDobl_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(QuadDobl_Complex_Ring,
                                QuadDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over quad double complex numbers.
