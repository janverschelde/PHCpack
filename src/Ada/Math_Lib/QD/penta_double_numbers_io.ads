with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;

package Penta_Double_Numbers_io is

-- DESCRIPTION :
--   This package defines input and output of penta double numbers.
--   The code is based on QD-2.3.9 (Y. Hida, X.S. Li, and D.H. Bailey).

  procedure read ( s : in string; ps : in out integer;
                   x : out penta_double; fail : out boolean );

  -- DESCRIPTION :
  --   Reads the string s starting at the position ps.

  -- ON ENTRY :
  --   s        contains an penta double number;
  --   ps       position in string s where to start reading.

  -- ON RETURN :
  --   ps       position in s of next valid entry (in case of a pair);
  --   x        penta double read from the string;
  --   fail     if format on string s is not correct.

  procedure read ( s : in string; x : out penta_double; fail : out boolean );

  -- DESCRIPTION :
  --    Reads a penta double from the string s, the result is stored in x.
  --    If all goes well, fail is false on return,
  --    otherwise fail is true on return.

  procedure to_digits ( x : in penta_double; precision : in natural32;
                        s : out string; expn : out integer32 );

  -- DESCRIPTION :
  --    Converts the penta double in x to decimal format,
  --    written to the string s, called by to_string below.
  --    See the specification of to_string for more information.

  procedure to_string ( x : in penta_double; precision,width : in natural32;
                        fixed,showpos,uppercase : in boolean;
                        fill : in character;
                        s : out string; endstring : out natural );

  -- DESCRIPTION :
  --   Writes the penta double x to the string in s.

  -- ON ENTRY :
  --   x          a penta double;
  --   s          a string large enough to hold at least d+8 characters,
  --              where d is the number of significant digits to write;
  --   precision  equals the number of significant digits to write;
  --   width      if large enough, the string s will be padded with
  --              characters equals to fill;
  --   fixed      if 1, then fixed format is used,
  --              otherwise, x is displayed in scientific format;
  --   showpos    if 1, then the + sign will be shown for positive numbers;
  --   uppercase  if 1, then E, otherwise e is used for the exponent;
  --   fill       character used for padding the string in case width
  --              is large enough. 
 
  -- ON RETURN :
  --   s          string that represents the penta double x;
  --   endstring  marks the end of the string: s = s(1..endstring).

  procedure write ( x : in penta_double; precision : in natural32 );

  -- DESCRIPTION :
  --    Writes the penta double to standard output, with given precision,
  --    and using default values for the to_string.

  procedure get ( x : in out penta_double );
  procedure get ( file : in file_type; x : in out penta_double );

  -- DESCRIPTION :
  --   Reads a string of characters from standard input or from file 
  --   and then reads the string into a penta double returned in x.

  procedure get ( x,y : in out penta_double );
  procedure get ( file : in file_type; x,y : in out penta_double );

  -- DESCRIPTION :
  --   Reads a string from file and then reads the string into two 
  --   penta doubles x and y.  This input routine is useful to read
  --   the real and imaginary part of a penta double complex number.

  procedure put ( x : in penta_double );
  procedure put ( x : in penta_double; precision : in natural32 );
  procedure put ( file : in file_type; x : in penta_double );
  procedure put ( file : in file_type; x : in penta_double;
                  precision : in natural32 );

  -- DESCRIPTION :
  --   Writes the penta double x to standard output or to file.
  --   By default, precision equals 80.

end Penta_Double_Numbers_io;
