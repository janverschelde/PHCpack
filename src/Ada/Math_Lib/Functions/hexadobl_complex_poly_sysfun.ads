with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_VecVecs;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Poly_Functions;
with HexaDobl_Complex_Poly_Systems;
with Generic_Poly_System_Functions;

package HexaDobl_Complex_Poly_SysFun is
  new Generic_Poly_System_Functions(HexaDobl_Complex_Ring,
                                    HexaDobl_Complex_Vectors,
                                    HexaDobl_Complex_VecVecs,
                                    HexaDobl_Complex_Polynomials,
                                    HexaDobl_Complex_Poly_Functions,
                                    HexaDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials over the
--   ring of hexa double complex numbers.
