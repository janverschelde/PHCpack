with TripDobl_Complex_Ring;
with TripDobl_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package TripDobl_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(TripDobl_Complex_Ring,
                                TripDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over triple double complex numbers.
