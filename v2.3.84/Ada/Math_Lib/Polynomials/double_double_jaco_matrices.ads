with Double_Double_Ring;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with Double_Double_Matrices;
with Double_Double_Polynomials;
with Double_Double_Poly_Functions;
with Double_Double_Poly_Systems;
with Double_Double_Poly_SysFun;
with Generic_Jacobian_Matrices;

package Double_Double_Jaco_Matrices is
  new Generic_Jacobian_Matrices(Double_Double_Ring,
                                Double_Double_Vectors,
                                Double_Double_VecVecs,
                                Double_Double_Matrices,
                                Double_Double_Polynomials,
                                Double_Double_Poly_Functions,
                                Double_Double_Poly_Systems,
                                Double_Double_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials with double double coefficients.
