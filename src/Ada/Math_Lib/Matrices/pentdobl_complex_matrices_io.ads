with PentDobl_Complex_Ring_io;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Matrices;
with Generic_Matrices_io;

package PentDobl_Complex_Matrices_io is 
  new Generic_Matrices_io(PentDobl_Complex_Ring_io,
                          PentDobl_Complex_Vectors,
                          PentDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of penta double complex numbers.
