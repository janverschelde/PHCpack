with Multprec_Complex_Ring;
with Multprec_Complex_Ring.FField;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laur_Functions;
with Multprec_Complex_Laur_Systems;
with Generic_Laur_System_Functions;

package Multprec_Complex_Laur_SysFun is
  new Generic_Laur_System_Functions(Multprec_Complex_Ring,
                                    Multprec_Complex_Ring.FField,
                                    Multprec_Complex_Vectors,
                                    Multprec_Complex_VecVecs,
                                    Multprec_Complex_Laurentials,
                                    Multprec_Complex_Laur_Functions,
                                    Multprec_Complex_Laur_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of Laurent polynomials
--   over the  multiprecision complex numbers.
