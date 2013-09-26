with Multprec_Natural_Ring_io;
with Multprec_Natural_Vectors;
with Multprec_Natural_Vectors_io;
with Multprec_Natural_VecVecs;
with Generic_VecVecs_io;

package Multprec_Natural_VecVecs_io is 
  new Generic_VecVecs_io(Multprec_Natural_Ring_io,
                         Multprec_Natural_Vectors,
                         Multprec_Natural_Vectors_io,
                         Multprec_Natural_VecVecs);

-- DESCRIPTION :
--   Defines i/o of vectors of vectors of multi-precision natural numbers.
