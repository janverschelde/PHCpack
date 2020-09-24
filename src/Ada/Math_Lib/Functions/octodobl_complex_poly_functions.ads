with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Polynomials;
with Generic_Polynomial_Functions;

package OctoDobl_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(OctoDobl_Complex_Ring,
                                   OctoDobl_Complex_Vectors,
                                   OctoDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for complex octo double numbers.
