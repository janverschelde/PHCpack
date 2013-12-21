with DoblDobl_Complex_Ring_io;
with DoblDobl_Complex_Vectors;
with Generic_Vectors_io;

package DoblDobl_Complex_Vectors_io is 
  new Generic_Vectors_io(DoblDobl_Complex_Ring_io,DoblDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of double double complex numbers.
