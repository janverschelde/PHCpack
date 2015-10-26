with Quad_Double_Ring;
with Quad_Double_Vectors;
with Quad_Double_VecVecs;
with Quad_Double_Polynomials;
with Quad_Double_Poly_Functions;
with Quad_Double_Poly_Systems;
with Generic_Poly_System_Functions;

package Quad_Double_Poly_SysFun is
  new Generic_Poly_System_Functions(Quad_Double_Ring,
                                    Quad_Double_Vectors,
                                    Quad_Double_VecVecs,
                                    Quad_Double_Polynomials,
                                    Quad_Double_Poly_Functions,
                                    Quad_Double_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials 
--   with real quad double coefficients.
