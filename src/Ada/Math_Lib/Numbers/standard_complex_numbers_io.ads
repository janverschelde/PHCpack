with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Standard_Complex_Numbers_io is

-- DESCRIPTION :
--   This package defines basic input/output for standard complex numbers.
--   In the basic format, the real and imaginary part of the complex
--   numbers are separated by spaces.

  procedure get ( c : in out Complex_Number );
  procedure get ( file : in file_type; c : in out Complex_Number );
  procedure get ( s : in string; c : in out Complex_Number;
                  last : out integer );

  -- DESCRIPTION :
  --   Reads a complex number c from standard output, from file,
  --   or from string s.

  procedure put ( c : in Complex_Number );
  procedure put ( file : in file_type; c : in Complex_Number );
  procedure put ( s : out string; c : in Complex_Number );

  -- DESCRIPTION :
  --   Writes the complex number c in its default format
  --   to standard output, to file, or to a string.

  procedure put ( c : in Complex_Number; fore,aft,exp : in natural32 );
  procedure put ( file : in file_type;
                  c : in Complex_Number; fore,aft,exp : in natural32 );
  procedure put ( s : out string;
                  c : in Complex_Number; aft,exp : in natural32 );

  -- DESCRIPTION :
  --   Writes a complex number to standard output, to file, or string,
  --   with the specified format.

  procedure put ( c : in Complex_Number; dp : in natural32 );
  procedure put ( file : in file_type;
                  c : in Complex_Number; dp : in natural32 );
  procedure put ( s : out string;
                  c : in Complex_Number; dp : in natural32 );

  -- DESCRIPTION :
  --   The complex number c is written with as many decimal places
  --   as the value for dp, to standard output, to file, or to string.

end Standard_Complex_Numbers_io;
