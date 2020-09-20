with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Polynomials;
with Generic_Polynomial_Systems;

package OctoDobl_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(OctoDobl_Complex_Ring,
                                 OctoDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the octo double complex numbers.
