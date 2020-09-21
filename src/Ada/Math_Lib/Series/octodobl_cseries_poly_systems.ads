with OctoDobl_Complex_Series_Ring;
with OctoDobl_CSeries_Polynomials;
with Generic_Polynomial_Systems;

package OctoDobl_CSeries_Poly_Systems is
  new Generic_Polynomial_Systems(OctoDobl_Complex_Series_Ring,
                                 OctoDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines systems of polynomials in several variables with coefficients
--   as series of octo double precision complex numbers.
