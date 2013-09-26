with Multprec_Floating_Ring_io;
with Multprec_Floating_Vectors;
with Multprec_Floating_Vectors_io;
with Multprec_Floating_VecVecs;
with Generic_VecVecs_io;

package Multprec_Floating_VecVecs_io is 
  new Generic_VecVecs_io(Multprec_Floating_Ring_io,
                         Multprec_Floating_Vectors,
                         Multprec_Floating_Vectors_io,
                         Multprec_Floating_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of multi-precision floating
--   numbers.
