with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Symbol_Table;                       use Symbol_Table;

package Line_Breaks is

  right : constant natural32 := 49;

  procedure Init_Line;

  -- DESCRIPTION :
  --   Initializes the internal column counter to zero.

  procedure Line ( file : in file_type; n : natural32 );

  -- DESCRIPTION :
  --   Takes a new line if n plus the current value of column
  --   exceeds the value of right, where n is the number of characters
  --   that will be written on the output file.

  function Length ( d,i : integer32; std : boolean; pw : power )
                  return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters needed to write one factor.

  -- ON ENTRY :
  --   d        degree of the factor;
  --   i        number of the variable;
  --   std      true when symbol table is not used;
  --   pw       symbol used for exponentiation.

end Line_Breaks;
