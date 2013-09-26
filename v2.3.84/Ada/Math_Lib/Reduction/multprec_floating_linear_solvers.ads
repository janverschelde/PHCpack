with Multprec_Floating_Ring;
with Multprec_Floating_Ring.FField;
with Multprec_Floating_Vectors;
with Multprec_Floating_Matrices;
with Generic_Floating_LInear_Solvers;

package Multprec_Floating_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Multprec_Floating_Ring,
                                      Multprec_Floating_Ring.FField,
                                      Multprec_Floating_Vectors,
                                      Multprec_Floating_Matrices);

-- DESCRIPTION :
--   Defines solvers for multi-precision floating-point linear systems.
