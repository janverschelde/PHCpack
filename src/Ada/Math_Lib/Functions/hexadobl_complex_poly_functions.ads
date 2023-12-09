with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_Polynomials;
with Generic_Polynomial_Functions;

package HexaDobl_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(HexaDobl_Complex_Ring,
                                   HexaDobl_Complex_Vectors,
                                   HexaDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for complex hexa double numbers.
