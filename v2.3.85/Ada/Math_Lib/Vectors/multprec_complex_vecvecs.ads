with Generic_VecVecs;
with Multprec_Complex_Ring;
with Multprec_Complex_Vectors;

package Multprec_Complex_VecVecs is 
  new Generic_VecVecs(Multprec_Complex_Ring,Multprec_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the multi-precision complex numbers.
