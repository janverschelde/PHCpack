with Multprec_Natural_Ring_io;
with Multprec_Natural_Vectors;
with Generic_Vectors_io;

package Multprec_Natural_Vectors_io is 
  new Generic_Vectors_io(Multprec_Natural_Ring_io,Multprec_Natural_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors over the ring of multi-precision
--   natural numbers.
