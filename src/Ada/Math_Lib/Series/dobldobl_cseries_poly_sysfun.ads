with DoblDobl_Complex_Series_Ring;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Functions;
with DoblDobl_CSeries_Poly_Systems;
with Generic_Poly_System_Functions;

package DoblDobl_CSeries_Poly_SysFun is
  new Generic_Poly_System_Functions(DoblDobl_Complex_Series_Ring,
                                    DoblDobl_Complex_Series_Vectors,
                                    DoblDobl_Complex_Series_VecVecs,
                                    DoblDobl_CSeries_Polynomials,
                                    DoblDobl_CSeries_Poly_Functions,
                                    DoblDobl_CSeries_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   double double precision complex numbers as coefficients.
