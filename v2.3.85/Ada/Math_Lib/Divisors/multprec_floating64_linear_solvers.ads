with Multprec_Floating64_Ring;
with Multprec_Floating64_Ring.FField;
with Multprec_Floating64_Vectors;
with Multprec_Floating64_Matrices;
with Generic_Floating_Linear_Solvers;

package Multprec_Floating64_Linear_Solvers is
  new Generic_Floating_Linear_Solvers(Multprec_Floating64_Ring,
                                      Multprec_Floating64_Ring.FField,
                                      Multprec_Floating64_Vectors,
                                      Multprec_Floating64_Matrices);

-- DESCRIPTION :
--   Defines solvers for multi-precision floating-point linear systems.
