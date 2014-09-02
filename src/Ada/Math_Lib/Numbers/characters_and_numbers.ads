with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Characters_and_Numbers is

-- DESCRIPTION :
--   This package offers some utilities for the input/output routines
--   of arbitrary precision numbers.

  function Integer_to_Character ( i : integer32 ) return character;
  function Character_to_Integer ( c : character ) return integer32;

  -- DESCRIPTION :
  --   The character corresponding to the integer code in i is returned,
  --   or vice versa.

  function Convert ( c : character ) return natural32;

  -- DESCRIPTION :
  --   Returns 10 if the character does not represent a number
  --   between 0 and 9, otherwise returns the corresponding number.

  function Convert ( s : string ) return natural32;
  function Convert ( s : string ) return natural64;

  -- DESCRIPTION :
  --   Converts the string into a 32-bit or 64-bit natural number,
  --   the string is supposed to contain a number in decimal format.

  function Convert ( s : string ) return integer32;
  function Convert ( s : string ) return integer64;

  -- DESCRIPTION :
  --   Converts the string into a 32-bit or 64-bit integer number,
  --   the string is supposed to contain a number in decimal format.

  function Convert ( s : string ) return double_float;

  -- DESCRIPTION :
  --   Converts the string representation of a natural number.
  --   The double float is there to deal with large numbers.

  function Convert_Decimal ( n : natural32 ) return character;

  -- DESCRIPTION :
  --   Returns the character representation of the number n in [0,9].

  function Convert_Hexadecimal ( c : character ) return natural32;

  -- DESCRIPTION :
  --   Returns the numerical representation of the character,
  --   if it is a hexadecimal symbol, otherwise 16 is returned.

  function Convert_Hexadecimal ( n : natural32 ) return character;

  -- DESCRIPTION :
  --   Returns the character representation of the number n in [0,15].

  function nConvert ( n : natural32 ) return string;
  function nConvert ( n : natural64 ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the natural number.

  function Convert ( i : integer32 ) return string;
  function Convert ( i : integer64 ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the integer number. 

  procedure Skip_Spaces ( file : in file_type; c : in out character );

  -- DESCRIPTION :
  --   Scans the file for the first character that is not a space.
  --   That character is returned as the parameter c on return.

  procedure Skip_Underscores ( file : in file_type; c : in out character );

  -- DESCRIPTION :
  --   Scans the file for the first character that is not an underscore.
  --   That character is returned as the parameter c on return.

end Characters_and_Numbers;
