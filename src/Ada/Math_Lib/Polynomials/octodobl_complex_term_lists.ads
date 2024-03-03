with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Polynomials;
with Generic_Lists_of_Terms;

package OctoDobl_Complex_Term_Lists is
  new Generic_Lists_of_Terms(OctoDobl_Complex_Ring,
                             OctoDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms for polynomial with complex coefficients 
--   in octo double precision.
