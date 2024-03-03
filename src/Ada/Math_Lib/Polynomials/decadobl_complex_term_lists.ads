with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Polynomials;
with Generic_Lists_of_Terms;

package DecaDobl_Complex_Term_Lists is
  new Generic_Lists_of_Terms(DecaDobl_Complex_Ring,
                             DecaDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms for polynomial with complex coefficients 
--   in deca double precision.
