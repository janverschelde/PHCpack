with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package DoblDobl_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(DoblDobl_Complex_Ring,
                                DoblDobl_Complex_Vectors,
                                DoblDobl_Complex_VecVecs,
                                DoblDobl_Complex_Matrices,
                                DoblDobl_Complex_Polynomials,
                                DoblDobl_Complex_Poly_Functions,
                                DoblDobl_Complex_Poly_Systems,
                                DoblDobl_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the double double complex numbers.
