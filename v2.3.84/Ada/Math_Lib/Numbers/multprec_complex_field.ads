with Multprec_Floating_Ring;
with Multprec_Floating_Ring.Field;
with Multprec_Complex_Ring;
with Multprec_Complex_Ring.Field;
with Generic_Complex_Field;

package Multprec_Complex_Field is
  new Generic_Complex_Field(Multprec_Floating_Ring,
                            Multprec_Floating_Ring.Field,
                            Multprec_Complex_Ring,
                            Multprec_Complex_Ring.Field);

-- DESCRIPTION :
--   Defines the complex field of multi-precision numbers.
