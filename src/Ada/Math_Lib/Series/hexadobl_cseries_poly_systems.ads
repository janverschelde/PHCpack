with HexaDobl_Complex_Series_Ring;
with HexaDobl_CSeries_Polynomials;
with Generic_Polynomial_Systems;

package HexaDobl_CSeries_Poly_Systems is
  new Generic_Polynomial_Systems(HexaDobl_Complex_Series_Ring,
                                 HexaDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines systems of polynomials in several variables with coefficients
--   as series of hexa double precision complex numbers.
