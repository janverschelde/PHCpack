with Multprec_Complex_Ring_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;
with Multprec_Complex_VecVecs;
with Multprec_Complex_VecLists;
with Generic_Lists_of_Vectors_io;

package Multprec_Complex_VecLists_io is
  new Generic_Lists_of_Vectors_io(Multprec_Complex_Ring_io,
                                  Multprec_Complex_Vectors,
                                  Multprec_Complex_Vectors_io,
                                  Multprec_Complex_VecVecs,
                                  Multprec_Complex_VecLists);

-- DESCRIPTION :
--   Defines input/output for lists of multiprecision complex vectors.
