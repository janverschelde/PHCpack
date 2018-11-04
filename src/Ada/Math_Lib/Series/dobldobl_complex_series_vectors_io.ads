with DoblDobl_Complex_Series_Ring_io;
with DoblDobl_Complex_Series_Vectors;
with Generic_Vectors_io;

package DoblDobl_Complex_Series_Vectors_io is 
  new Generic_Vectors_io(DoblDobl_Complex_Series_Ring_io,
                         DoblDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of truncated power series
--   with as coefficients complex numbers in double double precision.
