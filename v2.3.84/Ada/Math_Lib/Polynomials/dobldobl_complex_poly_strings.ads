with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;      use DoblDobl_Complex_Poly_Systems;

package DoblDobl_Complex_Poly_Strings is

-- DESCRIPTION :
--   This package converts polynomials in several variables and with
--   double double complex floating-point coefficients into strings,
--   and vice versa.

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

end DoblDobl_Complex_Poly_Strings;
