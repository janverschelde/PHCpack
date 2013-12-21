with Quad_Double_Ring;
with Quad_Double_Ring.FField;
with Quad_Double_Vectors;
with Quad_Double_Matrices;
with Generic_Floating_LInear_Solvers;

package Quad_Double_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Quad_Double_Ring,
                                      Quad_Double_Ring.FField,
                                      Quad_Double_Vectors,
                                      Quad_Double_Matrices);

-- DESCRIPTION :
--   Defines solvers for linear systems of quad doubles.
