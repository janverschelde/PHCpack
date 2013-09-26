with Multprec_Complex_Ring_io;
with Multprec_Complex_Vectors;
with Generic_Vectors_io;

package Multprec_Complex_Vectors_io is 
  new Generic_Vectors_io(Multprec_Complex_Ring_io,Multprec_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of multi-precision complex numbers.
