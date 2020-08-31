with Penta_Double_Ring;
with Penta_Double_Ring.FField;
with Penta_Double_Vectors;
with Penta_Double_Matrices;
with Generic_Floating_LInear_Solvers;

package Penta_Double_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Penta_Double_Ring,
                                      Penta_Double_Ring.FField,
                                      Penta_Double_Vectors,
                                      Penta_Double_Matrices);

-- DESCRIPTION :
--   Defines solvers for linear systems of penta doubles.
