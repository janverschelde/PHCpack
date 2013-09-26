with Multprec_Integer64_Ring_io;
with Multprec_Integer64_Vectors;
with Generic_Vectors_io;

package Multprec_Integer64_Vectors_io is 
  new Generic_Vectors_io(Multprec_Integer64_Ring_io,
                         Multprec_Integer64_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors over the ring of multi-precision
--   integer numbers, using 64-bit arithmetic.
