with OctoDobl_Complex_Series_Ring_io;
with OctoDobl_Complex_Series_Vectors;
with Generic_Vectors_io;

package OctoDobl_Complex_Series_Vectors_io is 
  new Generic_Vectors_io(OctoDobl_Complex_Series_Ring_io,
                         OctoDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of truncated power series
--   with as coefficients complex numbers in octo double precision.
