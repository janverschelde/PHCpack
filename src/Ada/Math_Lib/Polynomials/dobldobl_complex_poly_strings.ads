with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Symbol_Table;                       use Symbol_Table;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Term_Lists;        use DoblDobl_Complex_Term_Lists;
with DoblDobl_Complex_Poly_Systems;      use DoblDobl_Complex_Poly_Systems;

package DoblDobl_Complex_Poly_Strings is

-- DESCRIPTION :
--   This package converts polynomials in several variables and with
--   double double complex floating-point coefficients into strings,
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
  --   starting at s(k).  In case of a list, p_last points to
  --   the last term in the term list p.

  function Parse ( n : natural32; s : string ) return Poly;
  function Parse ( n : natural32; s : string ) return Term_List;

  -- DESCRIPTION :
  --   This function returns a polynomial in n variables,
  --   as a result of parsing the string s.

  function Parse ( n,m : natural32; s : string ) return Poly_Sys;
  function Parse ( n,m : natural32; s : string ) return Array_of_Term_Lists;

  -- DESCRIPTION :
  --   Returns a system of n polynomials in m variables.

  function Parse ( m : natural32; s : Array_of_Strings ) return Poly_Sys;
  function Parse ( m : natural32; s : Array_of_Strings )
                 return Array_of_Term_Lists;

  -- DESCRIPTION :
  --   Parses the strings in s into polynomials in m variables.

  function Size_Limit ( p : Poly ) return natural32;

  -- DESCRIPTION :
  --   Estimates the size of the string representation of the polynomial,
  --   providing a likely upper bound as needed for allocation in the
  --   interface programs, multiplying the number of terms with the 
  --   number of variables, the coefficient size and the symbol size.
  --   The size is limited by 2**32 - 1, the largest positive integer.

  function Write ( p : Poly ) return string;
  function Write ( p : Poly; s : Array_of_Symbols ) return string;

  -- DESCRIPTION :
  --   This function writes the polynomial to a string.
  --   Without s, the symbols in the symbol table represent the variables,
  --   otherwise, with s, the variables are written with the symbols in s.

  function Write ( p : Poly_Sys ) return string;
  function Write ( p : Poly_Sys; s : Array_of_Symbols ) return string;

  -- DESCRIPTION :
  --   This function writes the polynomial system to a string.
  --   Without s, the symbols in the symbol table represent the variables,
  --   otherwise, with s, the variables are written with the symbols in s.

  function Write ( p : Poly_Sys ) return Array_of_Strings;
  function Write ( p : Poly_Sys; s : Array_of_Symbols )
                 return Array_of_Strings;

  -- DESCRIPTION :
  --   Writes every polynomial in p to a separate string.
  --   Without s, the symbols in the symbol table represent the variables,
  --   otherwise, with s, the variables are written with the symbols in s.

end DoblDobl_Complex_Poly_Strings;
