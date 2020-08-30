with OctoDobl_Complex_Ring_io;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Matrices;
with Generic_Matrices_io;

package OctoDobl_Complex_Matrices_io is 
  new Generic_Matrices_io(OctoDobl_Complex_Ring_io,
                          OctoDobl_Complex_Vectors,
                          OctoDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of octo double complex numbers.
