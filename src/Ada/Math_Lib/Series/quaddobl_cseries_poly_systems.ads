with QuadDobl_Complex_Series_Ring;
with QuadDobl_CSeries_Polynomials;
with Generic_Polynomial_Systems;

package QuadDobl_CSeries_Poly_Systems is
  new Generic_Polynomial_Systems(QuadDobl_Complex_Series_Ring,
                                 QuadDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines systems of polynomials in several variables with coefficients
--   as series of quad double precision complex numbers.
