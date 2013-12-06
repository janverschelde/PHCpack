with Multprec_Floating64_Ring;
with Generic_Vectors;

package Multprec_Floating64_Vectors is
  new Generic_Vectors(Multprec_Floating64_Ring);

-- DESCRIPTION :
--   Defines vectors over the ring of multi-precision floating-point numbers,
--   for 64-bit arithmetic in fraction and exponent.
