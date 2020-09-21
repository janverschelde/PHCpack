with PentDobl_Complex_Series_Ring;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Complex_Series_VecVecs;
with PentDobl_CSeries_Polynomials;
with PentDobl_CSeries_Poly_Functions;
with PentDobl_CSeries_Poly_Systems;
with Generic_Poly_System_Functions;

package PentDobl_CSeries_Poly_SysFun is
  new Generic_Poly_System_Functions(PentDobl_Complex_Series_Ring,
                                    PentDobl_Complex_Series_Vectors,
                                    PentDobl_Complex_Series_VecVecs,
                                    PentDobl_CSeries_Polynomials,
                                    PentDobl_CSeries_Poly_Functions,
                                    PentDobl_CSeries_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   penta double precision complex numbers as coefficients.
