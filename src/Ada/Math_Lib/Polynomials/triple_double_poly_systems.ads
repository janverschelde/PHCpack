with Triple_Double_Ring;
with Triple_Double_Polynomials;
with Generic_Polynomial_Systems;

package Triple_Double_Poly_Systems is
  new Generic_Polynomial_Systems(Triple_Double_Ring,
                                 Triple_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with triple double coefficients.
