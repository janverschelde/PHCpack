with QuadDobl_Complex_Ring_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Generic_Matrices_io;

package QuadDobl_Complex_Matrices_io is 
  new Generic_Matrices_io(QuadDobl_Complex_Ring_io,
                          QuadDobl_Complex_Vectors,
                          QuadDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of quad double complex numbers.
