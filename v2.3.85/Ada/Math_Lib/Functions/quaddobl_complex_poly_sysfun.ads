with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Systems;
with Generic_Poly_System_Functions;

package QuadDobl_Complex_Poly_SysFun is
  new Generic_Poly_System_Functions(QuadDobl_Complex_Ring,
                                    QuadDobl_Complex_Vectors,
                                    QuadDobl_Complex_VecVecs,
                                    QuadDobl_Complex_Polynomials,
                                    QuadDobl_Complex_Poly_Functions,
                                    QuadDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials over the
--   ring of quad double complex numbers.
