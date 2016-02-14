with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with String_Splitters;                   use String_Splitters;

package Parse_Dimensions is

-- DESCRIPTION :
--   Given one polynomial as a string, or a list of strings,
--   the functions in this package parse the strings and return
--   the number of variables that occur in the polynomials.

  function Dim ( maxvar : natural32; strpol : string ) return natural32;

  -- DESCRIPTION :
  --   Given in maxvar the maximum number of variables that may occur
  --   in the string strpol, returns the number of variables in strpol.
  --   The symbol table is initialized with maxvar and the symbols
  --   that occur are still present after the function terminates.
  --   The string gets parsed as a Laurent polynomial with standard
  --   double precision complex coefficients.
  --   The parsed polynomial is available with the Get function
  --   and can be cleared with the Clear procedure.

  function Dim ( maxvar : natural32; strsys : Array_of_Strings )
               return natural32;

  -- DESCRIPTION :
  --   Given in maxvar the maximum number of variables that may occur
  --   in the strings strsys, returns the number of variables in strsys.
  --   The symbol table is initialized with maxvar and the symbols that
  --   occur in the Laurent polynomial systems are still present after
  --   the termination of the function.  The strings are parsed as Laurent
  --   polynomials with standard double precision coefficients.
  --   The parsed Laurent system is available with the Get function
  --   and can be cleared with the Clear procedure.

  function Get return Poly;

  -- DESCRIPTION :
  --   Returns the parsed polynomial.

  function Get return Link_to_Laur_Sys;

  -- DESCRIPTION :
  --   Returns the parsed system.

  procedure Clear;

  -- DESCRIPTION :
  --   Deallocates the space for the parsed polynomial and system.

end Parse_Dimensions;
