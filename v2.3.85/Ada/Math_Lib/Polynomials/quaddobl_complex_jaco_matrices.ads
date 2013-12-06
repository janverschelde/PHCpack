with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package QuadDobl_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(QuadDobl_Complex_Ring,
                                QuadDobl_Complex_Vectors,
                                QuadDobl_Complex_VecVecs,
                                QuadDobl_Complex_Matrices,
                                QuadDobl_Complex_Polynomials,
                                QuadDobl_Complex_Poly_Functions,
                                QuadDobl_Complex_Poly_Systems,
                                QuadDobl_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the quad double complex numbers.
