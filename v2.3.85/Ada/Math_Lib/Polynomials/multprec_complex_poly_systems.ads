with Multprec_Complex_Ring;
with Multprec_Complex_Polynomials;
with Generic_Polynomial_Systems;

package Multprec_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(Multprec_Complex_Ring,
                                 Multprec_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the multi-precision complex numbers.
