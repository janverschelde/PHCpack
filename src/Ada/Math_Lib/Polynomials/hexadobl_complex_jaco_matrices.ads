with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_VecVecs;
with HexaDobl_Complex_Matrices;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Poly_Functions;
with HexaDobl_Complex_Poly_Systems;
with HexaDobl_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package HexaDobl_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(HexaDobl_Complex_Ring,
                                HexaDobl_Complex_Vectors,
                                HexaDobl_Complex_VecVecs,
                                HexaDobl_Complex_Matrices,
                                HexaDobl_Complex_Polynomials,
                                HexaDobl_Complex_Poly_Functions,
                                HexaDobl_Complex_Poly_Systems,
                                HexaDobl_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the hexa double complex numbers.
