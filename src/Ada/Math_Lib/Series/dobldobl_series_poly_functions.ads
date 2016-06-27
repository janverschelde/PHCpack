with DoblDobl_Dense_Series_Ring;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Series_Polynomials;
with Generic_Polynomial_Functions;

package DoblDobl_Series_Poly_Functions is
  new Generic_Polynomial_Functions(DoblDobl_Dense_Series_Ring,
                                   DoblDobl_Dense_Series_Vectors,
                                   DoblDobl_Series_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with double double complex coefficients.
