with HexaDobl_Complex_Series_Ring_io;
with HexaDobl_Complex_Series_Vectors;
with Generic_Vectors_io;

package HexaDobl_Complex_Series_Vectors_io is 
  new Generic_Vectors_io(HexaDobl_Complex_Series_Ring_io,
                         HexaDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of truncated power series
--   with as coefficients complex numbers in hexa double precision.
