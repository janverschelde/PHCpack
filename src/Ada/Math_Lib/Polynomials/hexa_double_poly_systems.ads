with Hexa_Double_Ring;
with Hexa_Double_Polynomials;
with Generic_Polynomial_Systems;

package Hexa_Double_Poly_Systems is
  new Generic_Polynomial_Systems(Hexa_Double_Ring,
                                 Hexa_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with hexa double coefficients.
