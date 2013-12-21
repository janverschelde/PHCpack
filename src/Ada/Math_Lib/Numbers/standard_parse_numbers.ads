with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Standard_Parse_Numbers is

-- DESCRIPTION :
--    Numbers are parsed from file character by character.
--    The numbers may be given as rational numbers like 2/3, or 1.5,
--    or even using scientific notation for floating-point numbers.
--    Since hardware numbers are returned, there are limitations
--    on the length of the numbers.

  INFINITE_NUMBER : exception;
      -- occurs when a rational coefficient has a denominator = 0

  function Is_to_Skip ( ch : character ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the the charater is a space, carriage return,
  --   or line feed.

  procedure Skip_Spaces_and_CR
                  ( file : in file_type; ch : in out character );

  -- DESCRIPTION :
  --   As long as the character ch is a space, or carriage return, or
  --   line feed, a new character is read from file and returned in ch.
 
  -- ON RETURN :
  --   Either End_of_File(file) or ch is no space and no carriage return.

  procedure Skip_Spaces_and_CR
                  ( s : in string; p : in out integer );

  -- DESCRIPTION :
  --   As long as the character s(p) is a space or carriage return,
  --   the position p in the string is advanced.
 
  -- ON RETURN :
  --   Either p > s'last or s(p) is no space and no carriage return.

  procedure Parse ( file : in file_type;
                    char : in out character; i : out integer32;
                    ni : out natural32; sign : out character );
  procedure Parse_also_Brackets
                  ( file : in file_type;
                    char : in out character; i : out integer32;
                    ni : out natural32; sign : out character );

  -- DESCRIPTION :
  --   Characters are read from the input and a number is build up;
  --   the result is the integer number i.
  --   In the "also_Brackets", the number may be enclosed by round brackets.
  
  -- ON ENTRY :
  --   file         file type of a file that must be opened for input;
  --   char         first character to be analized.
 
  -- ON RETURN :
  --   char         first character that is not a digit;
  --   i1           digits read;
  --   ni           number of digits in i;
  --   sign         sign of the number.

  procedure Parse ( s : in string; p : in out integer;
                    i : out integer32; ni : out natural32;
                    sign : out character );
  procedure Parse_also_Brackets
                  ( s : in string; p : in out integer;
                    i : out integer32; ni : out natural32;
                    sign : out character );

  -- DESCRIPTION :
  --   The string is parsed for an integer number, starting at s(p).
  --   In the "also_Brackets", the number may be enclosed by round brackets.
  
  -- ON ENTRY :
  --   s            string of characters, positioned at p;
  --   p            current position in the string.
 
  -- ON RETURN :
  --   p            points at first nondigit character in s;
  --   i1           digits read;
  --   ni           number of digits in i;
  --   sign         sign of the number.

  procedure Parse ( file : in file_type;
                    char : in out character; i1,i2 : out integer32;
                    ni1,ni2 : out natural32; sign : out character );
  procedure Parse_also_Brackets
                  ( file : in file_type;
                    char : in out character; i1,i2 : out integer32;
                    ni1,ni2 : out natural32; sign : out character );

  -- DESCRIPTION :
  --   Characters are read from the input and a number is built;
  --   the result is the number : i1*10^ni2 + i2.
  --   In the "also_Brackets", the number may be enclosed by round brackets.
  
  -- ON ENTRY :
  --   file         file type of a file that must be opened for input;
  --   char         first character to be parsed.
 
  -- ON RETURN :
  --   char         first character that is not a digit;
  --   i1, i2       digits read;
  --   ni1, ni2     number of digits in i1 and i2;
  --   sign         sign of the number.

  procedure Parse ( s : in string; p : in out integer;
                    i1,i2 : out integer32; ni1,ni2 : out natural32;
                    sign : out character );
  procedure Parse_also_Brackets
                  ( s : in string; p : in out integer;
                    i1,i2 : out integer32; ni1,ni2 : out natural32;
                    sign : out character );

  -- DESCRIPTION :
  --   Parses the characters from the string and a number is built;
  --   the result is the number : i1*10^ni2 + i2.
  --   In the "also_Brackets", the number may be enclosed by round brackets.
  
  -- ON ENTRY :
  --   s            string of characters, s(p) is start of number;
  --   p            current position in the string.
 
  -- ON RETURN :
  --   p            s(p) is first character that is not a digit;
  --   i1, i2       digits read;
  --   ni1, ni2     number of digits in i1 and i2;
  --   sign         sign of the number.

  procedure Parse ( file : in file_type;
                    char : in out character; f : out double_float );

  -- DESCRIPTION :
  --   Parses a floating-point number, reading characters from file.

  procedure Parse ( s : in string; p : in out integer;
                    f : out double_float );

  -- DESCRIPTION :
  --   Parses a string into a floating-point number, starting at s(p).

  procedure Parse ( file : in file_type;
                    char : in out character; c : out Complex_Number );

  -- DESCRIPTION :
  --   A floating point number is read and converted into a complex number;
  --   the number may be the quotient of two floating point numbers.

  procedure Parse ( s : in string; p : in out integer;
                    c : out Complex_Number );

  -- DESCRIPTION :
  --   The string is scanned for a floating-point number that is then
  --   converted into a complex number, starting at the character s(p);
  --   The number may be the quotient of two floating point numbers.
 
end Standard_Parse_Numbers;
