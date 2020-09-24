with PentDobl_Complex_Ring;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Poly_Functions;
with PentDobl_Complex_Poly_Systems;
with Generic_Poly_System_Functions;

package PentDobl_Complex_Poly_SysFun is
  new Generic_Poly_System_Functions(PentDobl_Complex_Ring,
                                    PentDobl_Complex_Vectors,
                                    PentDobl_Complex_VecVecs,
                                    PentDobl_Complex_Polynomials,
                                    PentDobl_Complex_Poly_Functions,
                                    PentDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials over the
--   ring of penta double complex numbers.
