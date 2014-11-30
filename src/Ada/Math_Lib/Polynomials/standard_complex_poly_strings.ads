with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Standard_Complex_Poly_Strings is

-- DESCRIPTION :
--   This package converts polynomials in several variables and with
--   standard complex floating-point coefficients into strings, and
--   vice versa.

-- STRING MANIPULATORS :

  procedure Read_Exponent
              ( s : in string; k : in out integer; e : out natural32 );

  -- DESCRIPTION :
  --   Starts reading a natural number in the string at position k.

  -- ON ENTRY :
  --   s        string where s(k) is a digit;
  --   k        position to start the parsing.

  -- ON RETURN :
  --   k        position in s just after last digit read;
  --   e        value of the number read.

  function Delimiters ( n : natural32; s : string )
                      return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Searches the string s for n semicolons and returns the positions
  --   of those semicolons.

  function Concat_Symbol0 ( s : string; sb : Symbol ) return string;

  -- DESCRIPTION :
  --   Returns the result of concatenating the symbol to the string s.
  --   There is no '*' between the string and the symbol.

  function Concat_Symbol1 ( s : string; sb : Symbol ) return string;

  -- DESCRIPTION :
  --   Returns the result of concatenating the symbol to the string s.
  --   There is a '*' between the string and the symbol.

-- PARSING OPERATORS :

  procedure Parse ( s : in string; k : in out integer;
                    n : in natural32; p : in out Poly );

  -- DESCRIPTION :
  --   Parses the string s for a polynomial p in n variables,
  --   starting at s(k).

  function Parse ( n : natural32; s : string ) return Poly;

  -- DESCRIPTION :
  --   This function returns a polynomial in n variables,
  --   as a result of parsing the string s.

  function Parse ( n,m : natural32; s : string ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns a system of n polynomials in m variables.

  function Parse ( m : natural32; s : Array_of_Strings ) return Poly_Sys;

  -- DESCRIPTION :
  --   Parses the strings in s into polynomials in m variables.

  function Write ( p : Poly ) return string;

  -- DESCRIPTION :
  --   This function writes the polynomial to a string.

  function Write ( p : Poly_Sys ) return string;

  -- DESCRIPTION :
  --   This function writes the polynomial system to a string.

  function Write ( p : Poly_Sys ) return Array_of_Strings;

  -- DESCRIPTION :
  --   Writes each polynomial in p to a separate string.

end Standard_Complex_Poly_Strings;
