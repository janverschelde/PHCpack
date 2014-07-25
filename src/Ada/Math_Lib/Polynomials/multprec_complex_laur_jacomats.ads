with Multprec_Complex_Ring;
with Multprec_Complex_Ring.FField;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laur_Functions;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_SysFun;
with Generic_Laur_Jaco_Matrices;

package Multprec_Complex_Laur_JacoMats is
  new Generic_Laur_Jaco_Matrices(Multprec_Complex_Ring,
                                 Multprec_Complex_Ring.FField,
                                 Multprec_Complex_Vectors,
                                 Multprec_Complex_VecVecs,
                                 Multprec_Complex_Matrices,
                                 Multprec_Complex_Laurentials,
                                 Multprec_Complex_Laur_Functions,
                                 Multprec_Complex_Laur_Systems,
                                 Multprec_Complex_Laur_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of Laurent polynomials with multiprecision complex numbers.
