with Quad_Double_Ring;
with Quad_Double_Polynomials;
with Generic_Polynomial_Systems;

package Quad_Double_Poly_Systems is
  new Generic_Polynomial_Systems(Quad_Double_Ring,
                                 Quad_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with quad double coefficients.
