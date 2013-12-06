with Generic_NesVecs;
with Multprec_Complex_Ring;
with Multprec_Complex_Vectors;

package Multprec_Complex_NesVecs is 
  new Generic_NesVecs(Multprec_Complex_Ring,Multprec_Complex_Vectors);

-- DESCRIPTION :
--   Defines nested vectors over the ring of multi-precision complex numbers.
