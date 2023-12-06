with HexaDobl_Complex_Series_Ring;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_Complex_Series_VecVecs;
with HexaDobl_CSeries_Polynomials;
with HexaDobl_CSeries_Poly_Functions;
with HexaDobl_CSeries_Poly_Systems;
with Generic_Poly_System_Functions;

package HexaDobl_CSeries_Poly_SysFun is
  new Generic_Poly_System_Functions(HexaDobl_Complex_Series_Ring,
                                    HexaDobl_Complex_Series_Vectors,
                                    HexaDobl_Complex_Series_VecVecs,
                                    HexaDobl_CSeries_Polynomials,
                                    HexaDobl_CSeries_Poly_Functions,
                                    HexaDobl_CSeries_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   hexa double precision complex numbers as coefficients.
