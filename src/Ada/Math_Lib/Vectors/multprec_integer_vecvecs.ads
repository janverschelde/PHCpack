with Generic_VecVecs;
with Multprec_Integer_Ring;
with Multprec_Integer_Vectors;

package Multprec_Integer_VecVecs is 
  new Generic_VecVecs(Multprec_Integer_Ring,Multprec_Integer_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors of multi-precision integer numbers.
