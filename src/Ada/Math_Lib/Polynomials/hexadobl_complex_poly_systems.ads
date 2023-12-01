with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Polynomials;
with Generic_Polynomial_Systems;

package HexaDobl_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(HexaDobl_Complex_Ring,
                                 HexaDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the hexa double complex numbers.
