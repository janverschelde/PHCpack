with Generic_VecVecs;
with Multprec_Floating_Ring;
with Multprec_Floating_Vectors;

package Multprec_Floating_VecVecs is 
  new Generic_VecVecs(Multprec_Floating_Ring,Multprec_Floating_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors of multi-precision floating numbers.
