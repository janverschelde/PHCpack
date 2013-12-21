with Multprec_Complex_Ring;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Generic_Lists_of_Vectors;

package Multprec_Complex_VecLists is 
  new Generic_Lists_of_Vectors(Multprec_Complex_Ring,
                               Multprec_Complex_Vectors,
                               Multprec_Complex_VecVecs);

-- DESCRIPTION :
--   Defines lists of links to vectors of multiprecision complex numbers.
