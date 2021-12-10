with Hexa_Double_Ring;
with Hexa_Double_Ring.FField;
with Hexa_Double_Vectors;
with Hexa_Double_Matrices;
with Generic_Floating_LInear_Solvers;

package Hexa_Double_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Hexa_Double_Ring,
                                      Hexa_Double_Ring.FField,
                                      Hexa_Double_Vectors,
                                      Hexa_Double_Matrices);

-- DESCRIPTION :
--   Defines solvers for linear systems of hexa doubles.
