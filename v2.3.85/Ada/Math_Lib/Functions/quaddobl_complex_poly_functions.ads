with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with Generic_Polynomial_Functions;

package QuadDobl_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(QuadDobl_Complex_Ring,
                                   QuadDobl_Complex_Vectors,
                                   QuadDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for complex quad double numbers.
