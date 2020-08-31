with Triple_Double_Ring;
with Triple_Double_Ring.FField;
with Triple_Double_Vectors;
with Triple_Double_Matrices;
with Generic_Floating_LInear_Solvers;

package Triple_Double_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Triple_Double_Ring,
                                      Triple_Double_Ring.FField,
                                      Triple_Double_Vectors,
                                      Triple_Double_Matrices);

-- DESCRIPTION :
--   Defines solvers for linear systems of triple doubles.
