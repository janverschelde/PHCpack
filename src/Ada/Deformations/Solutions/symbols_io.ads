with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Symbol_Table;                       use Symbol_Table;

package Symbols_io is

-- DESCRIPTION :
--   Utilities to parse and write symbols, for the input and output
--   of solutions.

  procedure Skip_Symbol ( file : in file_type );

  -- DESCRIPTION :
  --   Skips all symbols until a `:' is encountered.

  function Read_Symbol ( file : file_type ) return Symbol;

  -- DESCRIPTION :
  --   Reads a symbol from file, skipping leading spaces, and returns it.

  function Get_Symbol ( file : file_type ) return natural32;

  -- DESCRIPTION :
  --   Reads a symbol from file and returns its number,
  --   returns 0 if the symbol does not occur in the table.

  function Get_Symbol ( file : file_type;
                        s : Array_of_Symbols ) return natural32;

  -- DESCRIPTION :
  --   Reads a symbol from file and returns its position in the array.

  procedure put_symbol ( file : in file_type; i : in natural32 );

  -- DESCRIPTION :
  --   Given the number of the symbol,
  --   the corresponding symbol will be written.

end Symbols_io;
