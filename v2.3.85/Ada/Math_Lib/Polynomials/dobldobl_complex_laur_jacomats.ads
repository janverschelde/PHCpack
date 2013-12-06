with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Ring.FField;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;
with Generic_Laur_Jaco_Matrices;

package DoblDobl_Complex_Laur_Jacomats is
  new Generic_Laur_Jaco_Matrices(DoblDobl_Complex_Ring,
                                 DoblDobl_Complex_Ring.FField,
                                 DoblDobl_Complex_Vectors,
                                 DoblDobl_Complex_VecVecs,
                                 DoblDobl_Complex_Matrices,
                                 DoblDobl_Complex_Laurentials,
                                 DoblDobl_Complex_Laur_Functions,
                                 DoblDobl_Complex_Laur_Systems,
                                 DoblDobl_Complex_Laur_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of Laurent polynomials over the double double complex numbers.
