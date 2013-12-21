with Generic_VecVecs;
with Multprec_Natural_Ring;
with Multprec_Natural_Vectors;

package Multprec_Natural_VecVecs is 
  new Generic_VecVecs(Multprec_Natural_Ring,Multprec_Natural_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors of multi-precision natural numbers.
