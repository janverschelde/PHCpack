with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Poly_Functions;
with OctoDobl_Complex_Poly_Systems;
with Generic_Poly_System_Functions;

package OctoDobl_Complex_Poly_SysFun is
  new Generic_Poly_System_Functions(OctoDobl_Complex_Ring,
                                    OctoDobl_Complex_Vectors,
                                    OctoDobl_Complex_VecVecs,
                                    OctoDobl_Complex_Polynomials,
                                    OctoDobl_Complex_Poly_Functions,
                                    OctoDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials over the
--   ring of octo double complex numbers.
