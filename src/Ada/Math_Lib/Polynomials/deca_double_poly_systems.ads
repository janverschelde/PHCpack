with Deca_Double_Ring;
with Deca_Double_Polynomials;
with Generic_Polynomial_Systems;

package Deca_Double_Poly_Systems is
  new Generic_Polynomial_Systems(Deca_Double_Ring,
                                 Deca_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with deca double coefficients.
