with Standard_Integer_Ring;
with Standard_Integer_Ring.DDomain;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Common_Divisors;
with Generic_Integer_LInear_Solvers;

package Standard_Integer_Linear_Solvers is
  new Generic_Integer_Linear_Solvers(Standard_Integer_Ring,
                                     Standard_Integer_Ring.DDomain,
                                     Standard_Integer_Vectors,
                                     Standard_Integer_Matrices,
                                     Standard_Common_Divisors);

-- DESCRIPTION :
--   Defines solvers for standard integer linear systems.
