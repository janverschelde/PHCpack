with PentDobl_Complex_Ring_io;
with PentDobl_Complex_Vectors;
with Generic_Vectors_io;

package PentDobl_Complex_Vectors_io is 
  new Generic_Vectors_io(PentDobl_Complex_Ring_io,PentDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of penta double complex numbers.
