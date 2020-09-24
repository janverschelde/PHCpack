with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Polynomials;
with Generic_Polynomial_Functions;

package DecaDobl_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(DecaDobl_Complex_Ring,
                                   DecaDobl_Complex_Vectors,
                                   DecaDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for complex deca double numbers.
