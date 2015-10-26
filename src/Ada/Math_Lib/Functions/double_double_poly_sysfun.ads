with Double_Double_Ring;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with Double_Double_Polynomials;
with Double_Double_Poly_Functions;
with Double_Double_Poly_Systems;
with Generic_Poly_System_Functions;

package Double_Double_Poly_SysFun is
  new Generic_Poly_System_Functions(Double_Double_Ring,
                                    Double_Double_Vectors,
                                    Double_Double_VecVecs,
                                    Double_Double_Polynomials,
                                    Double_Double_Poly_Functions,
                                    Double_Double_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials 
--   with real double double coefficients.
