with Standard_Complex_Ring;
with Standard_Complex_Ring.FField;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_Systems;
with Generic_Laur_System_Functions;

package Standard_Complex_Laur_SysFun is
  new Generic_Laur_System_Functions(Standard_Complex_Ring,
                                    Standard_Complex_Ring.FField,
                                    Standard_Complex_Vectors,
                                    Standard_Complex_VecVecs,
                                    Standard_Complex_Laurentials,
                                    Standard_Complex_Laur_Functions,
                                    Standard_Complex_Laur_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of Laurent polynomials over the
--   standard complex numbers.
