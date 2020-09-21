with DecaDobl_Complex_Series_Ring;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Complex_Series_VecVecs;
with DecaDobl_CSeries_Polynomials;
with DecaDobl_CSeries_Poly_Functions;
with DecaDobl_CSeries_Poly_Systems;
with Generic_Poly_System_Functions;

package DecaDobl_CSeries_Poly_SysFun is
  new Generic_Poly_System_Functions(DecaDobl_Complex_Series_Ring,
                                    DecaDobl_Complex_Series_Vectors,
                                    DecaDobl_Complex_Series_VecVecs,
                                    DecaDobl_CSeries_Polynomials,
                                    DecaDobl_CSeries_Poly_Functions,
                                    DecaDobl_CSeries_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   deca double precision complex numbers as coefficients.
