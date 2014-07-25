with Multprec_Complex_Ring;
with Multprec_Complex_Ring.FField;
with Multprec_Complex_Vectors;
with Multprec_Complex_Laurentials;
with Generic_Laur_Poly_Functions;

package Multprec_Complex_Laur_Functions is
  new Generic_Laur_Poly_Functions(Multprec_Complex_Ring,
                                  Multprec_Complex_Ring.FField,
                                  Multprec_Complex_Vectors,
                                  Multprec_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial functions for multiprecision complex numbers.
