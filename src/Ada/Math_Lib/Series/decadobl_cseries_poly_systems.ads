with DecaDobl_Complex_Series_Ring;
with DecaDobl_CSeries_Polynomials;
with Generic_Polynomial_Systems;

package DecaDobl_CSeries_Poly_Systems is
  new Generic_Polynomial_Systems(DecaDobl_Complex_Series_Ring,
                                 DecaDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines systems of polynomials in several variables with coefficients
--   as series of deca double precision complex numbers.
