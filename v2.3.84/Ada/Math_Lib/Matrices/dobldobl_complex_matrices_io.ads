with DoblDobl_Complex_Ring_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with Generic_Matrices_io;

package DoblDobl_Complex_Matrices_io is 
  new Generic_Matrices_io(DoblDobl_Complex_Ring_io,
                          DoblDobl_Complex_Vectors,
                          DoblDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of double double complex numbers.
