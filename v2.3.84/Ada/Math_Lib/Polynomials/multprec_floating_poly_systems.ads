with Multprec_Floating_Ring;
with Multprec_Floating_Polynomials;
with Generic_Polynomial_Systems;

package Multprec_Floating_Poly_Systems is
  new Generic_Polynomial_Systems(Multprec_Floating_Ring,
                                 Multprec_Floating_Polynomials);

-- DESCRIPTION :
--   Defines systems of multivariate polynomials
--   with multiprecision floating-point coefficients.
