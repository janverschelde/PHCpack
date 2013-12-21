with Multprec_Floating_Ring_io;
with Multprec_Floating_Vectors;
with Multprec_Floating_Matrices;
with Generic_Matrices_io;

package Multprec_Floating_Matrices_io is 
  new Generic_Matrices_io(Multprec_Floating_Ring_io,
                          Multprec_Floating_Vectors,
                          Multprec_Floating_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of multi-precision floating-point numbers.
