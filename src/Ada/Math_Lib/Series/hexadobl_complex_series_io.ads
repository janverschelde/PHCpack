with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with HexaDobl_Complex_Series;             use HexaDobl_Complex_Series;

package HexaDobl_Complex_Series_io is

-- DESCRIPTION :
--   The package encapsulates the i/o for the coefficient vectors.

  procedure get ( s : in out Series );
  procedure get ( file : in file_type; s : in out Series );

  -- DESCRIPTION :
  --   Reads as many complex numbers as the degree plus one.

  procedure get ( s : in out Link_to_Series; deg : in integer32 );
  procedure get ( file : in file_type;
                  s : in out Link_to_Series; deg : in integer32 );

  -- DESCRIPTION :
  --   Reads as many complex numbers as the degree plus one.

  procedure get ( s : in out Link_to_Series );
  procedure get ( file : in file_type; s : in out Link_to_Series );

  -- DESCRIPTION :
  --   Reads first the degree and then
  --   as many complex numbers as the degree plus one.

  procedure put ( s : in Series );
  procedure put ( file : in file_type; s : in Series );
  procedure put ( s : in Link_to_Series );
  procedure put ( file : in file_type; s : in Link_to_Series );

  -- DESCRIPTION :
  --   Writes the coefficient vector to file or standard output,
  --   one coefficient per line.

  procedure put ( s : in Series; dp : in natural32 );
  procedure put ( file : in file_type;
                  s : in Series; dp : in natural32 );
  procedure put ( s : in Link_to_Series; dp : in natural32 );
  procedure put ( file : in file_type;
                  s : in Link_to_Series; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the coefficient vector to file or standard output,
  --   one coefficient per line, displayed with dp decimal places.

end HexaDobl_Complex_Series_io;
