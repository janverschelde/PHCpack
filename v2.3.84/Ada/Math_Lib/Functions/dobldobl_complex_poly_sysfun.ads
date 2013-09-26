with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;
with Generic_Poly_System_Functions;

package DoblDobl_Complex_Poly_SysFun is
  new Generic_Poly_System_Functions(DoblDobl_Complex_Ring,
                                    DoblDobl_Complex_Vectors,
                                    DoblDobl_Complex_VecVecs,
                                    DoblDobl_Complex_Polynomials,
                                    DoblDobl_Complex_Poly_Functions,
                                    DoblDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials over the
--   ring of double double complex numbers.
