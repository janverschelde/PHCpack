with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Polynomials;
with Generic_Lists_of_Terms;

package DoblDobl_Complex_Term_Lists is
  new Generic_Lists_of_Terms(DoblDobl_Complex_Ring,
                             DoblDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms for polynomial with complex coefficients 
--   in double double precision.
