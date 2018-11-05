with DoblDobl_Complex_Series_Ring;
with DoblDobl_CSeries_Polynomials;
with Generic_Polynomial_Systems;

package DoblDobl_CSeries_Poly_Systems is
  new Generic_Polynomial_Systems(DoblDobl_Complex_Series_Ring,
                                 DoblDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines systems of polynomials in several variables with coefficients
--   as series of double double precision complex numbers.
