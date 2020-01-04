with DoblDobl_Complex_Ring_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;
with DoblDobl_Complex_VecMats;
with Generic_VecMats_io;

package DoblDobl_Complex_VecMats_io is 
  new Generic_VecMats_io(DoblDobl_Complex_Ring_io,
                         DoblDobl_Complex_Vectors,
                         DoblDobl_Complex_Matrices,
                         DoblDobl_Complex_Matrices_io,
                         DoblDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of complex numbers
--   in double double precision.
