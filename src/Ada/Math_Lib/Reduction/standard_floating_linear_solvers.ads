with Standard_Floating_Ring;
with Standard_Floating_Ring.FField;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Generic_Floating_LInear_Solvers;

package Standard_Floating_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Standard_Floating_Ring,
                                      Standard_Floating_Ring.FField,
                                      Standard_Floating_Vectors,
                                      Standard_Floating_Matrices);

-- DESCRIPTION :
--   Defines solvers for standard floating-point linear systems.
