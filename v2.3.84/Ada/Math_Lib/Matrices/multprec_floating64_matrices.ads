with Multprec_Floating64_Ring;            use Multprec_Floating64_Ring;
with Multprec_Floating64_Vectors;
with Generic_Matrices;

package Multprec_Floating64_Matrices is
  new Generic_Matrices(Multprec_Floating64_Ring,
                       Multprec_Floating64_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of multi-precision floating-point numbers,
--   numbers using 64-bit integer arithmetic in fraction and exponent.
