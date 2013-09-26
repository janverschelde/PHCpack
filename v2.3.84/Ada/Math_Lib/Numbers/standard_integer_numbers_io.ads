with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Standard_Integer_Numbers_io is

-- DESCRIPTION :
--   Definitions of the input/output routines for standard integers.

  procedure get ( i : in out integer32 );
  procedure get ( file : in file_type; i : in out integer32 );
  procedure get ( i : in out integer64 );
  procedure get ( file : in file_type; i : in out integer64 );

  -- DESCRIPTION :
  --   Reads an integer number from file.

  procedure put ( i : in integer32 );
  procedure put ( file : in file_type; i : in integer32 );
  procedure put ( i : in integer64 );
  procedure put ( file : in file_type; i : in integer64 );

  -- DESCRIPTION :
  --   Writes an integer number on file, with no spaces in front.

  procedure put ( i : in integer32; dp : in natural32 );
  procedure put ( file : in file_type; i : in integer32; dp : in natural32 );
  procedure put ( i : in integer64; dp : in natural32 );
  procedure put ( file : in file_type; i : in integer64; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the number on file, using at least dp decimal places.
  --   If the number needs less space in its display,
  --   then blanks are added in front of the number.

end Standard_Integer_Numbers_io;
