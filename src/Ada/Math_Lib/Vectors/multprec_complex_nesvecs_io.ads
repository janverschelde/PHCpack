with Multprec_Complex_Ring_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;
with Multprec_Complex_NesVecs;
with Generic_NesVecs_io;

package Multprec_Complex_NesVecs_io is 
  new Generic_NesVecs_io(Multprec_Complex_Ring_io,
                         Multprec_Complex_Vectors,
                         Multprec_Complex_Vectors_io,
                         Multprec_Complex_NesVecs);

-- DESCRIPTION :
--   Defines input/output of nested vectors multi-precision complex numbers.
