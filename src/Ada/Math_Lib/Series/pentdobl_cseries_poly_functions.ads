with PentDobl_Complex_Series_Ring;
with PentDobl_Complex_Series_Vectors;
with PentDobl_CSeries_Polynomials;
with Generic_Polynomial_Functions;

package PentDobl_CSeries_Poly_Functions is
  new Generic_Polynomial_Functions(PentDobl_Complex_Series_Ring,
                                   PentDobl_Complex_Series_Vectors,
                                   PentDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with penta double complex coefficients.
