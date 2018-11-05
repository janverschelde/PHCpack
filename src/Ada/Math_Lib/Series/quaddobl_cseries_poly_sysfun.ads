with QuadDobl_Complex_Series_Ring;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Functions;
with QuadDobl_CSeries_Poly_Systems;
with Generic_Poly_System_Functions;

package QuadDobl_CSeries_Poly_SysFun is
  new Generic_Poly_System_Functions(QuadDobl_Complex_Series_Ring,
                                    QuadDobl_Complex_Series_Vectors,
                                    QuadDobl_Complex_Series_VecVecs,
                                    QuadDobl_CSeries_Polynomials,
                                    QuadDobl_CSeries_Poly_Functions,
                                    QuadDobl_CSeries_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   quad double precision complex numbers as coefficients.
