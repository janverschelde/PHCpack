with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_Matrices;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Poly_Functions;
with OctoDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package OctoDobl_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(OctoDobl_Complex_Ring,
                                OctoDobl_Complex_Vectors,
                                OctoDobl_Complex_VecVecs,
                                OctoDobl_Complex_Matrices,
                                OctoDobl_Complex_Polynomials,
                                OctoDobl_Complex_Poly_Functions,
                                OctoDobl_Complex_Poly_Systems,
                                OctoDobl_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the octo double complex numbers.
