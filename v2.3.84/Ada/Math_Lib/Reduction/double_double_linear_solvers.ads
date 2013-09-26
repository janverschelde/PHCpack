with Double_Double_Ring;
with Double_Double_Ring.FField;
with Double_Double_Vectors;
with Double_Double_Matrices;
with Generic_Floating_LInear_Solvers;

package Double_Double_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Double_Double_Ring,
                                      Double_Double_Ring.FField,
                                      Double_Double_Vectors,
                                      Double_Double_Matrices);

-- DESCRIPTION :
--   Defines solvers for linear systems of double doubles.
