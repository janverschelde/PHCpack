with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Integer64_Numbers;         use Multprec_Integer64_Numbers;

package Multprec_Integer64_Numbers_io is

-- DESCRIPTION :
--   This package provides basic input/output routines for integer numbers
--   of arbitrary length.

  procedure get ( lc : in out character; i : in out Integer_Number );
  procedure get ( file : in file_type;
                  lc : in out character;  i : in out Integer_Number );

  -- DESCRIPTION :
  --   Reads a string of numbers and returns an integer number.
  --   The parameter lc is the leading character on entry.
  --   On return it is the last character read.

  procedure get ( i : in out Integer_Number );
  procedure get ( file : in file_type; i : in out Integer_Number );

  -- DESCRIPTION :
  --   Reads a string of numbers and returns an integer number.

  procedure get ( lc : in out character; i : in out Integer_Number;
                  nonnegative : out boolean );
  procedure get ( file : in file_type;
                  lc : in out character;  i : in out Integer_Number;
                  nonnegative : out boolean );

  -- DESCRIPTION :
  --   Reads a string of numbers and returns an integer number.
  --   The parameter lc is the leading character on entry.
  --   On return in lc it is the last character read.
  --   The parameter nonnegative is false if a minus sign was found.
  --   In this case, the sign of i will be negative also in the case
  --   where i equals zero.  This is useful for reading floats.

  procedure get ( s : in string; i : in out Integer_Number );

  -- DESCRIPTION :
  --   Reads an integer number from the string s.

  procedure put ( i : in Integer_Number );
  procedure put ( file : in file_type; i : in Integer_Number );

  -- DESCRIPTION :
  --   Writes the number on Standard Output or on file.

  procedure put ( s : out string; i : in Integer_Number );

  -- DESCRIPTION :
  --   Writes the number i to the string.
  -- REQUIRED :
  --   s'last = Decimal_Places(Unsigned(i)) + 1 if Negative(i),
  --          = Decimal_Places(Unsigned(i)) otherwise.

  procedure put ( i : in Integer_Number; dp : in natural32 );
  procedure put ( file : in file_type;
                  i : in Integer_Number; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the number on Standard Output or on file, using at least
  --   dp decimal places.  If the number needs less space in its display,
  --   then blanks are added in front of the number.

end Multprec_Integer64_Numbers_io;
