with Multprec_Complex_Ring;
with Multprec_Complex_Vectors;
with Multprec_Complex_Polynomials;
with Generic_Polynomial_Functions;

package Multprec_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(Multprec_Complex_Ring,
                                   Multprec_Complex_Vectors,
                                   Multprec_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for multi-precision complex numbers.
