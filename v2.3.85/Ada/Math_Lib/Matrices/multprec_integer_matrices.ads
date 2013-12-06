with Multprec_Integer_Ring;
with Multprec_Integer_Vectors;
with Generic_Matrices;

package Multprec_Integer_Matrices is
  new Generic_Matrices(Multprec_Integer_Ring,
                       Multprec_Integer_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of multi-precision integer numbers.
