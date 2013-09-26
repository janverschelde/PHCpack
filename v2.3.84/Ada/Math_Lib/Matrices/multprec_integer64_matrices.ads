with Multprec_Integer64_Ring;
with Multprec_Integer64_Vectors;
with Generic_Matrices;

package Multprec_Integer64_Matrices is
  new Generic_Matrices(Multprec_Integer64_Ring,
                       Multprec_Integer64_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of multi-precision integer numbers,
--   using 64-bit arithmetic.
