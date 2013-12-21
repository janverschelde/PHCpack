with Standard_Complex_Ring;
with Standard_Complex_Ring.FField;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;
with Generic_Laur_Jaco_Matrices;

package Standard_Complex_Laur_Jacomats is
  new Generic_Laur_Jaco_Matrices(Standard_Complex_Ring,
                                 Standard_Complex_Ring.FField,
                                 Standard_Complex_Vectors,
                                 Standard_Complex_VecVecs,
                                 Standard_Complex_Matrices,
                                 Standard_Complex_Laurentials,
                                 Standard_Complex_Laur_Functions,
                                 Standard_Complex_Laur_Systems,
                                 Standard_Complex_Laur_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of Laurent polynomials over the standard complex numbers.
