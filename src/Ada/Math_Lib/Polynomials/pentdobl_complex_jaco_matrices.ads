with PentDobl_Complex_Ring;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_Matrices;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Poly_Functions;
with PentDobl_Complex_Poly_Systems;
with PentDobl_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package PentDobl_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(PentDobl_Complex_Ring,
                                PentDobl_Complex_Vectors,
                                PentDobl_Complex_VecVecs,
                                PentDobl_Complex_Matrices,
                                PentDobl_Complex_Polynomials,
                                PentDobl_Complex_Poly_Functions,
                                PentDobl_Complex_Poly_Systems,
                                PentDobl_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the penta double complex numbers.
