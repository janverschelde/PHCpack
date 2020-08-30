with OctoDobl_Complex_Ring_io;
with OctoDobl_Complex_Vectors;
with Generic_Vectors_io;

package OctoDobl_Complex_Vectors_io is 
  new Generic_Vectors_io(OctoDobl_Complex_Ring_io,OctoDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of octo double complex numbers.
