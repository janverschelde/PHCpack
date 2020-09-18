with OctoDobl_Complex_Ring_io;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Matrices;
with OctoDobl_Complex_Matrices_io;
with OctoDobl_Complex_VecMats;
with Generic_VecMats_io;

package OctoDobl_Complex_VecMats_io is 
  new Generic_VecMats_io(OctoDobl_Complex_Ring_io,
                         OctoDobl_Complex_Vectors,
                         OctoDobl_Complex_Matrices,
                         OctoDobl_Complex_Matrices_io,
                         OctoDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of complex numbers
--   in octo double precision.
