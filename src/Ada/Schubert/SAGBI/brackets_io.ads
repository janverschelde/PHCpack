with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Brackets;                           use Brackets;

package Brackets_io is

-- DESCRIPTION :
--   This package provides input/output operations for brackets.

  procedure get ( b : in out Bracket );
  procedure get ( file : in file_type; b : in out Bracket );
  procedure get ( b : in out Bracket; sign : out integer32 );
  procedure get ( file : in file_type; b : in out Bracket;
                  sign : out integer32 );

  -- DESCRIPTION :
  --   Expects as many natural number as the length of b,
  --   either from Standard_Input or from file if specified.
  --   Optionally, the sign of the permutation to order the entries
  --   is returned.

  procedure put ( b : in Bracket );
  procedure put ( file : in file_type; b : in Bracket );

  -- DESCRIPTION :
  --   Writes down the entries in the bracket, either on
  --   Standard_Output or on file if specified.

end Brackets_io;
