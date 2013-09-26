with Quad_Double_Ring;
with Quad_Double_Vectors;
with Quad_Double_Polynomials;
with Generic_Polynomial_Functions;

package Quad_Double_Poly_Functions is
  new Generic_Polynomial_Functions(Quad_Double_Ring,
                                   Quad_Double_Vectors,
                                   Quad_Double_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for evaluation with quad doubles.
