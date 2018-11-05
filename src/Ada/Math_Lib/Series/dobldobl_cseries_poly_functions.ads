with DoblDobl_Complex_Series_Ring;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_CSeries_Polynomials;
with Generic_Polynomial_Functions;

package DoblDobl_CSeries_Poly_Functions is
  new Generic_Polynomial_Functions(DoblDobl_Complex_Series_Ring,
                                   DoblDobl_Complex_Series_Vectors,
                                   DoblDobl_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with double double complex coefficients.
