with Multprec_Integer64_Ring_io;
with Multprec_Integer64_Vectors;
with Multprec_Integer64_Matrices;
with Generic_Matrices_io;

package Multprec_Integer64_Matrices_io is 
  new Generic_Matrices_io(Multprec_Integer64_Ring_io,
                          Multprec_Integer64_Vectors,
                          Multprec_Integer64_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of multi-precision integer numbers.
