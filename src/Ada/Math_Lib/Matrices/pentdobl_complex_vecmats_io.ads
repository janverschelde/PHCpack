with PentDobl_Complex_Ring_io;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Matrices;
with PentDobl_Complex_Matrices_io;
with PentDobl_Complex_VecMats;
with Generic_VecMats_io;

package PentDobl_Complex_VecMats_io is 
  new Generic_VecMats_io(PentDobl_Complex_Ring_io,
                         PentDobl_Complex_Vectors,
                         PentDobl_Complex_Matrices,
                         PentDobl_Complex_Matrices_io,
                         PentDobl_Complex_VecMats);

-- DESCRIPTION :
--   Defines input/output of vectors of matrices of complex numbers
--   in penta double precision.
