with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Ring.FField;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_Systems;
with Generic_Laur_System_Functions;

package QuadDobl_Complex_Laur_SysFun is
  new Generic_Laur_System_Functions(QuadDobl_Complex_Ring,
                                    QuadDobl_Complex_Ring.FField,
                                    QuadDobl_Complex_Vectors,
                                    QuadDobl_Complex_VecVecs,
                                    QuadDobl_Complex_Laurentials,
                                    QuadDobl_Complex_Laur_Functions,
                                    QuadDobl_Complex_Laur_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of Laurent polynomial
--   over the quad double complex numbers.
