with Octo_Double_Ring;
with Octo_Double_Ring.FField;
with Octo_Double_Vectors;
with Octo_Double_Matrices;
with Generic_Floating_LInear_Solvers;

package Octo_Double_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Octo_Double_Ring,
                                      Octo_Double_Ring.FField,
                                      Octo_Double_Vectors,
                                      Octo_Double_Matrices);

-- DESCRIPTION :
--   Defines solvers for linear systems of octo doubles.
