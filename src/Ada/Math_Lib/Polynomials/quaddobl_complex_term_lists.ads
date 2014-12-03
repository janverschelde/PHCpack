with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Polynomials;
with Generic_Lists_of_Terms;

package QuadDobl_Complex_Term_Lists is
  new Generic_Lists_of_Terms(QuadDobl_Complex_Ring,
                             QuadDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms for polynomial with complex coefficients 
--   in quad double precision.
