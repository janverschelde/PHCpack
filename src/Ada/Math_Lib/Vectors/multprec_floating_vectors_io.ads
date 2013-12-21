with Multprec_Floating_Ring_io;
with Multprec_Floating_Vectors;
with Generic_Vectors_io;

package Multprec_Floating_Vectors_io is 
  new Generic_Vectors_io(Multprec_Floating_Ring_io,Multprec_Floating_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of multi-precision floating-point numbers.
