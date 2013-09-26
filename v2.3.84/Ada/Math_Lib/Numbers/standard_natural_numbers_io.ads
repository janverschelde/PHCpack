with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Standard_Natural_Numbers_io is

-- DESCRIPTION :
--   Definitions of the input/output routines for standard naturals.

  procedure get ( n : in out natural32 );
  procedure get ( file : in file_type; n : in out natural32 );
  procedure get ( n : in out natural64 );
  procedure get ( file : in file_type; n : in out natural64 );

  -- DESCRIPTION :
  --   Reads a natural number from file.

  procedure put ( n : in natural32 );
  procedure put ( file : in file_type; n : in natural32 );
  procedure put ( n : in natural64 );
  procedure put ( file : in file_type; n : in natural64 );

  -- DESCRIPTION :
  --   Writes a natural number to the file, with no spaces in front.

  procedure put ( n,dp : in natural32 );
  procedure put ( file : in file_type; n,dp : in natural32 );
  procedure put ( n : in natural64; dp : in natural32 );
  procedure put ( file : in file_type; n : in natural64; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the number on file, using at least dp decimal places.
  --   If the number needs less space in its display,
  --   then blanks are added in front of the number. 

end Standard_Natural_Numbers_io;
