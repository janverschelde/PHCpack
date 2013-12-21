with Double_Double_Ring;
with Double_Double_Vectors;
with Double_Double_Polynomials;
with Generic_Polynomial_Functions;

package Double_Double_Poly_Functions is
  new Generic_Polynomial_Functions(Double_Double_Ring,
                                   Double_Double_Vectors,
                                   Double_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for evaluation with double doubles.
