with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Ring.FField;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_Systems;
with Generic_Laur_System_Functions;

package DoblDobl_Complex_Laur_SysFun is
  new Generic_Laur_System_Functions(DoblDobl_Complex_Ring,
                                    DoblDobl_Complex_Ring.FField,
                                    DoblDobl_Complex_Vectors,
                                    DoblDobl_Complex_VecVecs,
                                    DoblDobl_Complex_Laurentials,
                                    DoblDobl_Complex_Laur_Functions,
                                    DoblDobl_Complex_Laur_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of Laurent polynomial
--   over the double double complex numbers.
