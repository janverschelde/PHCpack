with Deca_Double_Ring;
with Deca_Double_Ring.FField;
with Deca_Double_Vectors;
with Deca_Double_Matrices;
with Generic_Floating_LInear_Solvers;

package Deca_Double_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Deca_Double_Ring,
                                      Deca_Double_Ring.FField,
                                      Deca_Double_Vectors,
                                      Deca_Double_Matrices);

-- DESCRIPTION :
--   Defines solvers for linear systems of deca doubles.
