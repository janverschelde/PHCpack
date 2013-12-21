with Multprec_Complex_Ring_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;
with Generic_Matrices_io;

package Multprec_Complex_Matrices_io is 
  new Generic_Matrices_io(Multprec_Complex_Ring_io,
                          Multprec_Complex_Vectors,
                          Multprec_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of multi-precision complex numbers.
