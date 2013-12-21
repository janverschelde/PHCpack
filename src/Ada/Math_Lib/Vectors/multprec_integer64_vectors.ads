with Multprec_Integer64_Ring;
with Generic_Vectors;

package Multprec_Integer64_Vectors is 
  new Generic_Vectors(Multprec_Integer64_Ring);

-- DESCRIPTION :
--   Defines vectors over the ring of multi-precision integer numbers,
--   using 64-bit arithmetic.
