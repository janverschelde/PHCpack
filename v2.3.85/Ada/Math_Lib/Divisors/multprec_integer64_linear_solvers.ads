with Multprec_Integer64_Ring;
with Multprec_Integer64_Ring.DDomain;
with Multprec_Integer64_Vectors;
with Multprec_Integer64_Matrices;
with Multprec64_Common_Divisors;
with Generic_Integer_LInear_Solvers;

package Multprec_Integer64_Linear_Solvers is
  new Generic_Integer_Linear_Solvers(Multprec_Integer64_Ring,
                                     Multprec_Integer64_Ring.DDomain,
                                     Multprec_Integer64_Vectors,
                                     Multprec_Integer64_Matrices,
                                     Multprec64_Common_Divisors);

-- DESCRIPTION :
--   Defines solvers for multiprecision integer linear systems,
--   using 64-bit arithmetic for the multiprecision numbers.
