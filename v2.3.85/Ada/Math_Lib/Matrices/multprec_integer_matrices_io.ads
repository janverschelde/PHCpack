with Multprec_Integer_Ring_io;
with Multprec_Integer_Vectors;
with Multprec_Integer_Matrices;
with Generic_Matrices_io;

package Multprec_Integer_Matrices_io is 
  new Generic_Matrices_io(Multprec_Integer_Ring_io,
                          Multprec_Integer_Vectors,
                          Multprec_Integer_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of multi-precision integer numbers.
