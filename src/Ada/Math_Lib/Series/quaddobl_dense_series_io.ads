with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with QuadDobl_Dense_Series;               use QuadDobl_Dense_Series;

package QuadDobl_Dense_Series_io is

-- DESCRIPTION :
--   The package encapsulates the i/o for the coefficient vectors,
--   of truncated power series in quad double precision.

  procedure get ( s : in out Series );
  procedure get ( file : in file_type; s : in out Series );

  -- DESCRIPTION :
  --   Prompts for the degree (an integer number),
  --   followed by as many complex numbers as the degree plus one.

  procedure put ( s : in Series );
  procedure put ( file : in file_type; s : in Series );

  -- DESCRIPTION :
  --   Writes the coefficient vector to file or standard output,
  --   one coefficient per line.

  procedure put ( s : in Series; dp : in natural32 );
  procedure put ( file : in file_type;
                  s : in Series; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the coefficient vector to file or standard output,
  --   one coefficient per line, displayed with dp decimal places.

end QuadDobl_Dense_Series_io;
