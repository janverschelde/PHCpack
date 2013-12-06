with text_io;                            use text_io;
with Symbol_Table;                       use Symbol_Table;

package Symbol_Table_io is

-- DESCRIPTION :
--   This package provides i/o operations for symbols.

  procedure get ( sb : in out symbol );
  procedure get ( sb : in out symbol; delimiter : in character );
  procedure get ( ch : in out character; sb : in out symbol );
  procedure get ( ch : in out character; sb : in out symbol;
                  delimiter : in character );
  procedure get ( file : in file_type; sb : in out symbol );
  procedure get ( file : in file_type; sb : in out symbol;
                  delimiter : in character );
  procedure get ( file : in file_type;
                  ch : in out character; sb : in out symbol );
  procedure get ( file : in file_type;
                  ch : in out character; sb : in out symbol;
                  delimiter : in character );

  -- DESCRIPTION :
  --   Reads a symbol from standard input or from file.
  --   Skips blanks before and after the symbol.
  --   When a character is given as parameter, it can already
  --   contain the first character of the symbol.
  --   Reading stops when either the length of the symbol is reached,
  --   or when a blank or a delimiter is encountered.

  procedure put ( sb : in symbol );
  procedure put ( file : in file_type; sb : in symbol );

  -- DESCRIPTION :
  --   Writes a symbol on standard output or on file.
  --   Blanks after the symbol are not written.

  procedure Write;
  procedure Write ( file : in file_type );

  -- DESCRIPTION :
  --   Writes all symbols in the table on standard output or on file.
  --   Each symbol is preceded by one space.

  procedure Parse_Symbol ( file : in file_type; ch : in out character;
                           sb : out Symbol );

  -- DESCRIPTION :
  --   Parses the input file for a symbol that starts with ch.
  --   Raises INVALID_SYMBOL if the symbol is not valid.

  procedure Parse_Symbol ( s : in string; p : in out integer;
                           sb : out Symbol );

  -- DESCRIPTION :
  --   Parses the string for a symbol, starting at position p.
  --   Raises INVALID_SYMBOL if the symbol is not valid.

end Symbol_Table_io;
