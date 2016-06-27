with DoblDobl_Dense_Series_Ring;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_VecVecs;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Functions;
with DoblDobl_Series_Poly_Systems;
with Generic_Poly_System_Functions;

package DoblDobl_Series_Poly_SysFun is
  new Generic_Poly_System_Functions(DoblDobl_Dense_Series_Ring,
                                    DoblDobl_Dense_Series_Vectors,
                                    DoblDobl_Dense_Series_VecVecs,
                                    DoblDobl_Series_Polynomials,
                                    DoblDobl_Series_Poly_Functions,
                                    DoblDobl_Series_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   double double complex numbers as coefficients.
