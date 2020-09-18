with DecaDobl_Complex_Ring_io;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Matrices;
with DecaDobl_Complex_Matrices_io;
with DecaDobl_Complex_VecMats;
with Generic_VecMats_io;

package DecaDobl_Complex_VecMats_io is 
  new Generic_VecMats_io(DecaDobl_Complex_Ring_io,
                         DecaDobl_Complex_Vectors,
                         DecaDobl_Complex_Matrices,
                         DecaDobl_Complex_Matrices_io,
                         DecaDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of complex numbers
--   in deca double precision.
