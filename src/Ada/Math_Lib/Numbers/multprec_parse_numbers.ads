with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;

package Multprec_Parse_Numbers is

-- DESCRIPTION :
--    Numbers are parsed from file character by character.
--    The numbers may be given as rational numbers like 2/3, or 1.5,
--    or even using scientific notation for floating-point numbers.
--    Since multiprecision numbers are returned, there is no limitation
--    on the length of the numbers.

  INFINITE_NUMBER : exception;
      -- occurs when a rational coefficient has a denominator = 0

  procedure Parse ( file : in file_type; char : in out character;
                    i : out Integer_Number;
                    n : out natural32; sign : out character );
  procedure Parse_also_Brackets
                  ( file : in file_type; char : in out character;
                    i : out Integer_Number;
                    n : out natural32; sign : out character );

  -- DESCRIPTION :
  --   Characters are read from the input and a number is parsed.
  --   In the "also_Brackets", the number may be enclosed by round Brackets.
  
  -- ON ENTRY :
  --   file         file type of a file that must be opened for input;
  --   char         first character to be analized.
 
  -- ON RETURN :
  --   char         first character that is not a digit;
  --   i            integer number read;
  --   n            number of digits in i;
  --   sign         sign of the number.

  procedure Parse ( s : in string; p : in out integer;
                    i : out Integer_Number; n : out natural32;
                    sign : out character );
  procedure Parse_also_Brackets
                  ( s : in string; p : in out integer;
                    i : out Integer_Number; n : out natural32;
                    sign : out character );

  -- DESCRIPTION :
  --   The string is parsed for an integer number, starting at s(p).
  --   In the "also_Brackets", the number may be enclosed by round brackets.
  
  -- ON ENTRY :
  --   s            string of characters, positioned at p;
  --   p            current position in the string.
 
  -- ON RETURN :
  --   p            points at first nondigit character in s;
  --   i            integer number read;
  --   n            number of digits in i;
  --   sign         sign of the number.

  procedure Parse ( file : in file_type; char : in out character;
                    f : out Floating_Number );

  -- DESCRIPTION :
  --   A multiprecision floating-point number is parsed from file.

  procedure Parse ( s : in string; p : in out integer;
                    f : out Floating_Number );

  -- DESCRIPTION :
  --   Parses the string, starting at position p for a multiprecision
  --   floating-point number.

  procedure Parse ( file : in file_type; size : in natural32;
                    char : in out character; c : out Complex_Number );

  -- DESCRIPTION :
  --   A floating point number is read and converted into a complex number;
  --   the number may be a fraction of two floating-point numbers.
  --   The size determines the precision of evaluating the fraction.
 
  procedure Parse ( s : in string; size : in natural32; p : in out integer;
                    c : out Complex_Number );

  -- DESCRIPTION :
  --   Parses the string starting at position p for a multiprecision
  --   floating-point number that is then converted into a complex number.
  --   The number may be a fraction of two floating-point numbers and
  --   the size determines the precision of evaluation the fraction.
 
end Multprec_Parse_Numbers;
