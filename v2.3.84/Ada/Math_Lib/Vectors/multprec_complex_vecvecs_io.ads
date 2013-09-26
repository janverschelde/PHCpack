with Multprec_Complex_Ring_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;
with Multprec_Complex_VecVecs;
with Generic_VecVecs_io;

package Multprec_Complex_VecVecs_io is 
  new Generic_VecVecs_io(Multprec_Complex_Ring_io,
                         Multprec_Complex_Vectors,
                         Multprec_Complex_Vectors_io,
                         Multprec_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of multprec complex numbers.
