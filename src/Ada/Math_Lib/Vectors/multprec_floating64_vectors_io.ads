with Multprec_Floating64_Ring_io;
with Multprec_Floating64_Vectors;
with Generic_Vectors_io;

package Multprec_Floating64_Vectors_io is 
  new Generic_Vectors_io(Multprec_Floating64_Ring_io,
                         Multprec_Floating64_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of multi-precision floating-point numbers.
