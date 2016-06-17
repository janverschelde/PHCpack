with text_io;                             use text_io;
with Standard_Dense_Series;               use Standard_Dense_Series;

package Standard_Dense_Series_io is

-- DESCRIPTION :
--   The package encapsulates the i/o for the coefficient vectors.

  procedure get ( s : out Series );
  procedure get ( file : in file_type; s : out Series );

  -- DESCRIPTION :
  --   Prompts for the order (an integer number),
  --   followed by as many complex numbers as the order plus one.

  procedure put ( s : in Series );
  procedure put ( file : in file_type; s : in Series );

  -- DESCRIPTION :
  --   Writes the coefficient vector to file or standard output,
  --   one coefficient per line.

end Standard_Dense_Series_io;
