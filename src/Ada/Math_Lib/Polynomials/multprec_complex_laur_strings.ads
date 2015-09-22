with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Complex_Laurentials;       use Multprec_Complex_Laurentials;
with Multprec_Complex_Laur_Systems;      use Multprec_Complex_Laur_Systems;

package Multprec_Complex_Laur_Strings is

-- DESCRIPTION :
--   This package converts Laurent polynomials in several variables 
--   with multiprecision complex floating-point coefficients into strings,
--   and vice versa.

  procedure Parse ( s : in string; k : in out integer;
                    n,size : in natural32; p : in out Poly );

  -- DESCRIPTION :
  --   Parses the string s for a polynomial p in n variables,
  --   starting at s(k) using working precision set at size
  --   to evaluate fractions.

  function Parse ( n,size : natural32; s : string ) return Poly;

  -- DESCRIPTION :
  --   This function returns a polynomial in n variables,
  --   as a result of parsing the string s, using working precision
   --  set at size to evaluate fractions.

  function Parse ( n,m,size : natural32; s : string ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns a system of n polynomials in m variables,
  --   using working precision set at size to evaluate fractions.

  function Parse ( m,size : natural32; s : Array_of_Strings ) return Laur_Sys;

  -- DESCRIPTION :
  --   Parses the strings in s into polynomials in m variables,
  --   using working precision set at size to evaluate fractions.

  function Size_Limit ( p : Poly ) return natural32;

  -- DESCRIPTION :
  --   Estimates the size of the string representation of the polynomial,
  --   providing a likely upper bound as needed for allocation in the
  --   interface programs, multiplying the number of terms with the 
  --   number of variables, the coefficient size and the symbol size.
  --   The size is limited by 2**32 - 1, the largest positive integer.

  function Write ( p : Poly ) return string;

  -- DESCRIPTION :
  --   This function writes the polynomial to a string.

  function Write ( p : Laur_Sys ) return string;

  -- DESCRIPTION :
  --   This function writes the polynomial system to a string.

  function Write ( p : Laur_Sys ) return Array_of_Strings;

  -- DESCRIPTION :
  --   Writes every polynomial in p to a separate string.

end Multprec_Complex_Laur_Strings;
