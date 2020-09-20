with Penta_Double_Ring;
with Penta_Double_Polynomials;
with Generic_Polynomial_Systems;

package Penta_Double_Poly_Systems is
  new Generic_Polynomial_Systems(Penta_Double_Ring,
                                 Penta_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with penta double coefficients.
