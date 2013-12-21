with Quad_Double_Ring;
with Quad_Double_Vectors;
with Quad_Double_VecVecs;
with Quad_Double_Matrices;
with Quad_Double_Polynomials;
with Quad_Double_Poly_Functions;
with Quad_Double_Poly_Systems;
with Quad_Double_Poly_SysFun;
with Generic_Jacobian_Matrices;

package Quad_Double_Jaco_Matrices is
  new Generic_Jacobian_Matrices(Quad_Double_Ring,
                                Quad_Double_Vectors,
                                Quad_Double_VecVecs,
                                Quad_Double_Matrices,
                                Quad_Double_Polynomials,
                                Quad_Double_Poly_Functions,
                                Quad_Double_Poly_Systems,
                                Quad_Double_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials with quad double coefficients.
