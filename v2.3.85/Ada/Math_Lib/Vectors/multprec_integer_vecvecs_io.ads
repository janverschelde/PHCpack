with Multprec_Integer_Ring_io;
with Multprec_Integer_Vectors;
with Multprec_Integer_Vectors_io;
with Multprec_Integer_VecVecs;
with Generic_VecVecs_io;

package Multprec_Integer_VecVecs_io is 
  new Generic_VecVecs_io(Multprec_Integer_Ring_io,
                         Multprec_Integer_Vectors,
                         Multprec_Integer_Vectors_io,
                         Multprec_Integer_VecVecs);

-- DESCRIPTION :
--   Defines i/o of vectors of vectors of multi-precision integer numbers.
