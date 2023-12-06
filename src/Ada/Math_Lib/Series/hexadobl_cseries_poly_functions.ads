with HexaDobl_Complex_Series_Ring;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_CSeries_Polynomials;
with Generic_Polynomial_Functions;

package HexaDobl_CSeries_Poly_Functions is
  new Generic_Polynomial_Functions(HexaDobl_Complex_Series_Ring,
                                   HexaDobl_Complex_Series_Vectors,
                                   HexaDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with hexa double complex coefficients.
