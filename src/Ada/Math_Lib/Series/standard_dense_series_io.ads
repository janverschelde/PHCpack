with text_io;                             use text_io;
with Standard_Dense_Series;               use Standard_Dense_Series;

package Standard_Dense_Series_io is

-- DESCRIPTION :
--   The package encapsulates the i/o for the coefficient vectors.

  procedure put ( s : in Series );
  procedure put ( file : in file_type; s : in Series );

  -- DESCRIPTION :
  --   Writes the coefficient vector to file or standard output,
  --   one coefficient per line.

end Standard_Dense_Series_io;
