with PentDobl_Complex_Ring;
with PentDobl_Complex_Polynomials;
with Generic_Lists_of_Terms;

package PentDobl_Complex_Term_Lists is
  new Generic_Lists_of_Terms(PentDobl_Complex_Ring,
                             PentDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms for polynomial with complex coefficients 
--   in penta double precision.
