with Multprec_Floating64_Ring_io;
with Multprec_Floating64_Vectors;
with Multprec_Floating64_Matrices;
with Generic_Matrices_io;

package Multprec_Floating64_Matrices_io is 
  new Generic_Matrices_io(Multprec_Floating64_Ring_io,
                          Multprec_Floating64_Vectors,
                          Multprec_Floating64_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of multi-precision floating-point numbers.
