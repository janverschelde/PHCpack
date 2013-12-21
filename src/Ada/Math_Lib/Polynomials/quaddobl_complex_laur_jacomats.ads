with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Ring.FField;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_SysFun;
with Generic_Laur_Jaco_Matrices;

package QuadDobl_Complex_Laur_Jacomats is
  new Generic_Laur_Jaco_Matrices(QuadDobl_Complex_Ring,
                                 QuadDobl_Complex_Ring.FField,
                                 QuadDobl_Complex_Vectors,
                                 QuadDobl_Complex_VecVecs,
                                 QuadDobl_Complex_Matrices,
                                 QuadDobl_Complex_Laurentials,
                                 QuadDobl_Complex_Laur_Functions,
                                 QuadDobl_Complex_Laur_Systems,
                                 QuadDobl_Complex_Laur_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of Laurent polynomials over the quad double complex numbers.
