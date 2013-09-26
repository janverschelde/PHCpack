with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Natural_Coefficients;      use Multprec_Natural_Coefficients;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;

package Multprec_Natural_Numbers_io is

-- DESCRIPTION :
--   This package provides basic input/output routines for natural numbers
--   of arbitrary length.  The input is restricted to numbers whose size
--   does not exceed the length of one line.
--   To enhance readability of long numbers, underscores may be used to
--   separate blocks of numbers.

  procedure get ( lc : in out character; n : in out Natural_Number );
  procedure get ( file : in file_type;
                  lc : in out character;  n : in out Natural_Number );
  procedure get ( n : in out Natural_Number );
  procedure get ( file : in file_type; n : in out Natural_Number );

  -- DESCRIPTION :
  --   Reads a string of numbers and returns a natural number.
  --   The parameter lc is the leading character on entry.
  --   On return it is the last character that has been read.

  procedure get ( s : in string; n : in out Natural_Number );

  -- DESCRIPTION :
  --   Creates a new natural number from the string s.

  procedure put ( n : in Natural_Number );
  procedure put ( file : in file_type; n : in Natural_Number );

  -- DESCRIPTION :
  --   Writes the number on Standard Output or on file.

  procedure put ( s : out string; n : in Natural_Number );

  -- DESCRIPTION :
  --   Writes the string representation of n.
  -- REQUIRED : s'last = Decimal_Places(n);

  procedure put ( n : in Array_of_Naturals );
  procedure put ( file : in file_type; n : in Array_of_Naturals );

  -- DESCRIPTION :
  --   Writes the array as a natural number with all leading zeros.

  procedure put ( s : out string; n : in Array_of_Naturals );

  -- DESCRIPTION :
  --   Writes the string representation of the coefficients in n to s.
  -- REQUIRED :
  --   s'last = n'length*Multprec_Natural_Numbers.Exponent.

  procedure put ( n : in Natural_Number; dp : in natural32 );
  procedure put ( file : in file_type;
                  n : in Natural_Number; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the number on Standard Output or on file, using at least
  --   dp decimal places.  If the number needs less space in its display,
  --   then blanks are added in front of the number.        

end Multprec_Natural_Numbers_io;
