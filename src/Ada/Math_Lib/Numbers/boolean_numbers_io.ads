with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Boolean_Numbers_io is

-- DESCRIPTION :
--   Definitions of the input/output routines for booleans,
--   represented by 0 for false and 1 for true.

  procedure get ( b : in out boolean );
  procedure get ( file : in file_type; b : in out boolean );

  -- DESCRIPTION :
  --   Reads an integer number from standard input or from file,
  --   returns in b false if the read integer is zero,
  --   returns in b true if the read integer is one.

  procedure put ( b : in boolean );
  procedure put ( file : in file_type; b : in boolean );

  -- DESCRIPTION :
  --   Writes 0 if b is false, 1 if b is true, without spaces in front,
  --   to standard output or to file.

  procedure put ( b : in boolean; dp : in natural32 );
  procedure put ( file : in file_type; b : in boolean; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the number on file, using at least dp decimal places.
  --   If the number needs less space in its display,
  --   then blanks are added in front of the number.

end Boolean_Numbers_io;
