with Multprec_Natural_Ring;
with Multprec_Natural_Vectors;
with Generic_Matrices;

package Multprec_Natural_Matrices is
  new Generic_Matrices(Multprec_Natural_Ring,
                       Multprec_Natural_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of multi-precision natural numbers.
