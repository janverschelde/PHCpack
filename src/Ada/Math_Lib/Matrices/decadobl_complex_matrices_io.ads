with DecaDobl_Complex_Ring_io;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Matrices;
with Generic_Matrices_io;

package DecaDobl_Complex_Matrices_io is 
  new Generic_Matrices_io(DecaDobl_Complex_Ring_io,
                          DecaDobl_Complex_Vectors,
                          DecaDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of deca double complex numbers.
