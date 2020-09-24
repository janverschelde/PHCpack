with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Poly_Functions;
with DecaDobl_Complex_Poly_Systems;
with Generic_Poly_System_Functions;

package DecaDobl_Complex_Poly_SysFun is
  new Generic_Poly_System_Functions(DecaDobl_Complex_Ring,
                                    DecaDobl_Complex_Vectors,
                                    DecaDobl_Complex_VecVecs,
                                    DecaDobl_Complex_Polynomials,
                                    DecaDobl_Complex_Poly_Functions,
                                    DecaDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials over the
--   ring of deca double complex numbers.
