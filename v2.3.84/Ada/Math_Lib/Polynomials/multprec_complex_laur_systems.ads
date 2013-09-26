with Multprec_Complex_Ring;
with Multprec_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package Multprec_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(Multprec_Complex_Ring,
                                Multprec_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over multi-precision complex numbers.
