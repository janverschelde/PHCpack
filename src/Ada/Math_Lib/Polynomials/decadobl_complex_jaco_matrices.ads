with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_Matrices;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Poly_Functions;
with DecaDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package DecaDobl_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(DecaDobl_Complex_Ring,
                                DecaDobl_Complex_Vectors,
                                DecaDobl_Complex_VecVecs,
                                DecaDobl_Complex_Matrices,
                                DecaDobl_Complex_Polynomials,
                                DecaDobl_Complex_Poly_Functions,
                                DecaDobl_Complex_Poly_Systems,
                                DecaDobl_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the deca double complex numbers.
