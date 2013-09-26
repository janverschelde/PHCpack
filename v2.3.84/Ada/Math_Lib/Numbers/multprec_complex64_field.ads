with Multprec_Floating64_Ring;
with Multprec_Floating64_Ring.Field;
with Multprec_Complex64_Ring;
with Multprec_Complex64_Ring.Field;
with Generic_Complex_Field;

package Multprec_Complex64_Field is
  new Generic_Complex_Field(Multprec_Floating64_Ring,
                            Multprec_Floating64_Ring.Field,
                            Multprec_Complex64_Ring,
                            Multprec_Complex64_Ring.Field);

-- DESCRIPTION :
--   Defines the complex field of multi-precision numbers.
