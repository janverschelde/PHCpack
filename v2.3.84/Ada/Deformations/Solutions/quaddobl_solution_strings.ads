with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;
with Symbol_Table;                       use Symbol_Table;

package QuadDobl_Solution_Strings is

-- DESCRIPTION :
--   This package provides facilities to convert solutions into strings
--   and to parse strings into solutions.
--   A string to represent a solution consists of an introduction
--   with the t and m information, the solution vector, and then
--   the diagnostics.
--   The operations are supposed to be idempotent: 
--   (1) solution -> write -> parse -> write -> same solution and
--   (2) string -> parse -> write -> parse -> same string.

-- PART I: compute lengths and write the strings

  function Length_Number ( f : quad_double ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters in the string
  --   representation of the floating-point number.

  function Length_Number ( c : Complex_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters in the string returned
  --   by Write_Number(c).

  function Write_Number ( c : Complex_Number ) return string;

  -- DESCRIPTION :
  --   Writes the complex number to string.

  function Length_Intro ( s : Solution ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters in the introduction
  --   of the string representation of the solution s.

  function Write_Intro ( s : Solution ) return string;

  -- DESCRIPTION :
  --   Writes the introduction to the solution s to string.

  function Length_Symbol ( i : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters need to represent the i-th symbol.

  function Write_Symbol ( i : natural32 ) return string;

  -- DESCRIPTION :
  --   Writes the symbol for the i-th variable to a string.

  function Length_Vector ( v : Vector ) return natural32;
  function Length_Vector ( s : Solution ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters in the string representation
  --   of the solution vector.

  function Write_Vector ( v : Vector ) return string;
  function Write_Vector ( s : Solution ) return string;

  -- DESCRIPTION :
  --   Writes the solution vector of s to a string.

  function Length_Diagnostics return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters in the string returned by
  --   the function Write_Diagnostics.

  function Write_Diagnostics ( err,rco,res : quad_double ) return string;
  function Write_Diagnostics ( s : Solution ) return string;

  -- DESCRIPTION :
  --   Writes the diagnostics line of a solution to a string.

  function Length ( s : Solution ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters in the string representation
  --   of the solution s.

  function Write ( s : Solution ) return string;

  -- DESCRIPTION :
  --   Writes a solution to a string in the same format
  --   as it is written to standard output or to file.

-- PART II: parse strings into solutions
--   Note that all exceptions are suppressed, i.e: handled silently with
--   default zero values on return.
--   If the symbol table is empty, then it is initialized with
--   the symbols in the strings.  If the number of symbols is less
--   than the dimension n, then additional symbols are added.

  procedure Parse_Intro
               ( s : in string; k : in out integer;
                 t : out Complex_Number; m : out integer32;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Parses the string for the value for the continuation parameter t
  --   and the multiplicity m.  If an exception has occurred, because
  --   the string s was not in the right format, then fail is true,
  --   otherwise fail is false.

  -- ON ENTRY :
  --   s         string containing solutions;
  --   k         index in s to the start of the introduction to a solution.
  
  -- ON RETURN :
  --   k         updated current position in the string,
  --             at the end of reading the introduction;
  --   t         is the parsed value for the continuation parameter;
  --   m         is the parsed value for multiplicity field;
  --   fail      if true, then an exception occurred and the data
  --             returned may be wrong, otherwise, it should be okay.
 
  procedure Parse_Symbol ( s : in String; k : in out integer;
                           sb : out Symbol; fail : out boolean );

  -- DESCRIPTION :
  --   Parses the string for a symbol, starting at the current
  --   position k in the string s, skipping over newline symbols
  --   and spaces.  On return in sb is the symbol between s(k)
  --   and the first occurrence of the colon in the string s.

  -- ON ENTRY :
  --   s         a string containing a solution vector;
  --   k         current position in the string, at start of symbol
  --             we expect newline symbol or space;

  -- ON RETURN :
  --   k         new current position in s where s(k) equals ':';
  --   sb        symbol read from the string;
  --   fail      true if some exception occurred, false otherwise.
 
  procedure Parse_Vector
               ( s : in string; k : in out integer; n : in natural32;
                 v : out QuadDobl_Complex_Vectors.Vector;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Parses the string for a solution vector of n components,
  --   starting at position k in the string.
  --   If fail is not false on return, then v contains the
  --   values of an n-dimensional solution vector.
  --   The current position k is updated.

  procedure Parse_Diagnostics
               ( s : in string; k : in out integer;
                 err,rco,res : out quad_double; fail : out boolean );

  -- DESCRIPTION :
  --   Parses the string for the solution diagnostics,
  --   starting at the current position k in the string.
  --   Updates the position k.
  --   If an exception occurs because s is in the wrong format,
  --   then fail will be true, otherwise fail is false.

  procedure Parse ( s : in string; k : in out integer; n : natural32;
                    sol : out Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Parses the string into a solution of n variables.

  -- ON ENTRY :
  --   s         a string with solutions;
  --   k         current position in the string,
  --             at the start of the introduction to a solution;
  --   n         number of variables in the solution.

  -- ON RETURN :
  --   k         updated current position in the string;
  --   sol       solution parsed from the string;
  --   fail      if true, then the format of the string was wrong
  --             or some exception occurred, otherwise, all is well.

end QuadDobl_Solution_Strings;
