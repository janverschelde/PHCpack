with Multprec_Complex_Ring_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;
with Multprec_Complex_Matrices_io;
with Multprec_Complex_VecMats;
with Generic_VecMats_io;

package Multprec_Complex_VecMats_io is 
  new Generic_VecMats_io(Multprec_Complex_Ring_io,
                         Multprec_Complex_Vectors,
                         Multprec_Complex_Matrices,
                         Multprec_Complex_Matrices_io,
                         Multprec_Complex_VecMats);

-- DESCRIPTION :
--   Input/output of vectors of matrices of multiprecision complex numbers.
