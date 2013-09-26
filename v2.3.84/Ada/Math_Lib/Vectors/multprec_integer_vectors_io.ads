with Multprec_Integer_Ring_io;
with Multprec_Integer_Vectors;
with Generic_Vectors_io;

package Multprec_Integer_Vectors_io is 
  new Generic_Vectors_io(Multprec_Integer_Ring_io,Multprec_Integer_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors over the ring of multi-precision
--   integer numbers.
