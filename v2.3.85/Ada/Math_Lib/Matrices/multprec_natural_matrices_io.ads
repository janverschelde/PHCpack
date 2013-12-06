with Multprec_Natural_Ring_io;
with Multprec_Natural_Vectors;
with Multprec_Natural_Matrices;
with Generic_Matrices_io;

package Multprec_Natural_Matrices_io is 
  new Generic_Matrices_io(Multprec_Natural_Ring_io,
                          Multprec_Natural_Vectors,
                          Multprec_Natural_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of multi-precision natural numbers.
