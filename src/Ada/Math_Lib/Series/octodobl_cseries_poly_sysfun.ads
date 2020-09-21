with OctoDobl_Complex_Series_Ring;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Complex_Series_VecVecs;
with OctoDobl_CSeries_Polynomials;
with OctoDobl_CSeries_Poly_Functions;
with OctoDobl_CSeries_Poly_Systems;
with Generic_Poly_System_Functions;

package OctoDobl_CSeries_Poly_SysFun is
  new Generic_Poly_System_Functions(OctoDobl_Complex_Series_Ring,
                                    OctoDobl_Complex_Series_Vectors,
                                    OctoDobl_Complex_Series_VecVecs,
                                    OctoDobl_CSeries_Polynomials,
                                    OctoDobl_CSeries_Poly_Functions,
                                    OctoDobl_CSeries_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   octo double precision complex numbers as coefficients.
