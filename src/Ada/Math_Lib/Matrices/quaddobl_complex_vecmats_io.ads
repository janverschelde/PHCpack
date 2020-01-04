with QuadDobl_Complex_Ring_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;
with QuadDobl_Complex_VecMats;
with Generic_VecMats_io;

package QuadDobl_Complex_VecMats_io is 
  new Generic_VecMats_io(QuadDobl_Complex_Ring_io,
                         QuadDobl_Complex_Vectors,
                         QuadDobl_Complex_Matrices,
                         QuadDobl_Complex_Matrices_io,
                         QuadDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of complex numbers
--   in quad double precision.
