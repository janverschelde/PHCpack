with Standard_Integer64_Ring;
with Standard_Integer64_Ring.DDomain;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;
with Standard64_Common_Divisors;
with Generic_Integer_LInear_Solvers;

package Standard_Integer64_Linear_Solvers is
  new Generic_Integer_Linear_Solvers(Standard_Integer64_Ring,
                                     Standard_Integer64_Ring.DDomain,
                                     Standard_Integer64_Vectors,
                                     Standard_Integer64_Matrices,
                                     Standard64_Common_Divisors);

-- DESCRIPTION :
--   Defines solvers for standard integer linear systems,
--   using 64-bit integer arithmetic.
