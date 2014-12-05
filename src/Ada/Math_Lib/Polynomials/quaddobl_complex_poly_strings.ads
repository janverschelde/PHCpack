with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with QuadDobl_Complex_Polynomials;       use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Term_Lists;        use QuadDobl_Complex_Term_Lists;
with QuadDobl_Complex_Poly_Systems;      use QuadDobl_Complex_Poly_Systems;

package QuadDobl_Complex_Poly_Strings is

-- DESCRIPTION :
--   This package converts polynomials in several variables and with
--   quad double complex floating-point coefficients into strings,
--   and vice versa.  The parse operations of strings into lists of terms
--   are limited: no support for powers of polynomials and long bracketed
--   expressions, because the goal of parsing into term lists is to avoid
--   the sorting and polynomial operations that are needed with the
--   parsing of general polynomial expressions.

  procedure Parse ( s : in string; k : in out integer;
                    n : in natural32; p : in out Poly );
  procedure Parse ( s : in string; k : in out integer;
                    n : in natural32; p,p_last : in out Term_List );

  -- DESCRIPTION :
  --   Parses the string s for a polynomial p in n variables,
  --   starting at s(k).

  function Parse ( n : natural32; s : string ) return Poly;
  function Parse ( n : natural32; s : string ) return Term_List;

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
  --   Writes every polynomial in p to a separate string.

end QuadDobl_Complex_Poly_Strings;
