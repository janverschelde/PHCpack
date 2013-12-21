with Multprec_Floating_Ring;
with Multprec_Floating_Vectors;
with Multprec_Floating_Polynomials;
with Generic_Polynomial_Functions;

package Multprec_Floating_Poly_Functions is
  new Generic_Polynomial_Functions(Multprec_Floating_Ring,
                                   Multprec_Floating_Vectors,
                                   Multprec_Floating_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for multi-precision floating-point numbers.
