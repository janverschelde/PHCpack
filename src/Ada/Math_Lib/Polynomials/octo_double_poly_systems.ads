with Octo_Double_Ring;
with Octo_Double_Polynomials;
with Generic_Polynomial_Systems;

package Octo_Double_Poly_Systems is
  new Generic_Polynomial_Systems(Octo_Double_Ring,
                                 Octo_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with octo double coefficients.
