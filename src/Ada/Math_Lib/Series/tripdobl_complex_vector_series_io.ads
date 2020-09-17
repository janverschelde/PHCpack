with text_io;                           use text_io;
with TripDobl_Complex_Vector_Series;    use TripDobl_Complex_Vector_Series;

package TripDobl_Complex_Vector_Series_io is

-- DESCRPTION :
--   Provides very basic output for series with vector coefficients.

  procedure put ( v : in Vector );
  procedure put ( file : in file_type; v : in Vector );

  -- DESCRIPTION :
  --   Writes the coefficients of the vector series v to standard output,
  --   or to file.  The coefficients are separated by a new line.

end TripDobl_Complex_Vector_Series_io;
