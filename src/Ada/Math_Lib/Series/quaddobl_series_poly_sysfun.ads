with QuadDobl_Dense_Series_Ring;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_VecVecs;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Functions;
with QuadDobl_Series_Poly_Systems;
with Generic_Poly_System_Functions;

package QuadDobl_Series_Poly_SysFun is
  new Generic_Poly_System_Functions(QuadDobl_Dense_Series_Ring,
                                    QuadDobl_Dense_Series_Vectors,
                                    QuadDobl_Dense_Series_VecVecs,
                                    QuadDobl_Series_Polynomials,
                                    QuadDobl_Series_Poly_Functions,
                                    QuadDobl_Series_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   quad double complex numbers as coefficients.
