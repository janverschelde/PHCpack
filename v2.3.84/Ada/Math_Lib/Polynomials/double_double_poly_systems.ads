with Double_Double_Ring;
with Double_Double_Polynomials;
with Generic_Polynomial_Systems;

package Double_Double_Poly_Systems is
  new Generic_Polynomial_Systems(Double_Double_Ring,
                                 Double_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with double double coefficients.
