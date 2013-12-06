with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with Generic_Polynomial_Functions;

package DoblDobl_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(DoblDobl_Complex_Ring,
                                   DoblDobl_Complex_Vectors,
                                   DoblDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for complex double double numbers.
