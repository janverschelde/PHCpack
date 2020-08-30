with DecaDobl_Complex_Ring_io;
with DecaDobl_Complex_Vectors;
with Generic_Vectors_io;

package DecaDobl_Complex_Vectors_io is 
  new Generic_Vectors_io(DecaDobl_Complex_Ring_io,DecaDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of deca double complex numbers.
