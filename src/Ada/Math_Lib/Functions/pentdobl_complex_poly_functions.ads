with PentDobl_Complex_Ring;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Polynomials;
with Generic_Polynomial_Functions;

package PentDobl_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(PentDobl_Complex_Ring,
                                   PentDobl_Complex_Vectors,
                                   PentDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for complex penta double numbers.
