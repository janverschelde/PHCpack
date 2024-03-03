with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Polynomials;
with Generic_Lists_of_Terms;

package HexaDobl_Complex_Term_Lists is
  new Generic_Lists_of_Terms(HexaDobl_Complex_Ring,
                             HexaDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms for polynomial with complex coefficients 
--   in hexa double precision.
