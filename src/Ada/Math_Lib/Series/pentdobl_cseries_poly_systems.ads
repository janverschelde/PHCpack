with PentDobl_Complex_Series_Ring;
with PentDobl_CSeries_Polynomials;
with Generic_Polynomial_Systems;

package PentDobl_CSeries_Poly_Systems is
  new Generic_Polynomial_Systems(PentDobl_Complex_Series_Ring,
                                 PentDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines systems of polynomials in several variables with coefficients
--   as series of penta double precision complex numbers.
