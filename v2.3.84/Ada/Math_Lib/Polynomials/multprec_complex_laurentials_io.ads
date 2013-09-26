with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Symbol_Table;                       use Symbol_Table;
with Multprec_Complex_Laurentials;       use Multprec_Complex_Laurentials;

package Multprec_Complex_Laurentials_io is

-- DESCRIPTION :
--   This package provides input/output routines for Laurent polynomials,
--   polynomials with exponents that may be negative and coefficients of
--   arbitrary precision.

-- SETTING THE WORKING PRECISION :

  procedure Set_Working_Precision ( size : in natural32 );

  -- DESCRIPTION :
  --   Sets the working precision to evaluate fractions in the get
  --   procedures below.

-- THE INPUT OPERATIONS :

  procedure get ( n : out natural32; p : out Poly );
  procedure get ( file : in file_type; n : out natural32; p : out Poly );
  procedure get ( p : out Poly );
  procedure get ( file : in file_type; p : out Poly );

  -- DESCRIPTION :
  --   Reads a multivariate polynomial from standard input or from file.

  -- ON ENTRY :
  --   file     optional file, must be opened for input.

  -- ON RETURN :
  --   n        optional number of variables, if provided then n will
  --            first be read and used to initialize the symbol table,
  --            otherwise, if n is not provided, then the symbol table
  --            must contain already n symbols;
  --   p        Laurent polynomial in n variables.

-- OUTPUT OPERATIONS :

  procedure put ( file : in file_type; d,i : in integer32;
                  std : in boolean; pow : in Power );
  procedure put ( file : in file_type; d : in integer32;
                  sb : in Symbol; pow : in Power );

  -- DESCRIPTION :
  --   Writes one factor x^d to file, where x is the given symbol,

  -- ON ENTRY :
  --   file     must be opened for output;
  --   d        degree of the factor;
  --   i        number of the variable in the symbol table
  --            if the symbol is not given and not standard;
  --   std      indicates whether the standard x symbols are used
  --            used instead of the symbols stored in the symbol table;
  --   sb       symbol to be used for the variable;
  --   pow      either "^" or "**".

  procedure put ( file : in file_type; d : in Degrees;
                  std : in boolean; pow : in Power );
  procedure put ( file : in file_type; d : in Degrees;
                  s : in Array_of_Symbols; pow : in Power );

  -- DESCRIPTION :
  --   Only those factors with nonzero degree d are written to file,
  --   separated by a multiplication '*'.

  procedure put ( t : in Term );
  procedure put ( t : in Term; s : in Array_of_Symbols );
  procedure put ( file : in file_type; t : in Term );
  procedure put ( file : in file_type; t : in Term;
                  s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Writes the term t to standard output or to file.

  -- ON ENTRY :
  --   file     optional file, must be opened for output;
  --   t        term of Laurent polynomial in several variables;
  --   s        optional list of symbols for the variables.

  procedure put ( p : in Poly );
  procedure put ( p : in Poly; s : in Array_of_Symbols );
  procedure put ( file : in file_type; p : in Poly );
  procedure put ( file : in file_type; p : in Poly;
                  s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Writes a polynomial to standard output or to file.

  -- ON ENTRY :
  --   file     optional file, must be opened for output;
  --   p        Laurent polynomial in several variables;
  --   s        optional array of symbols for the variables.

  procedure put_line ( p : in Poly );
  procedure put_line ( p : in Poly; s : in Array_of_Symbols );
  procedure put_line ( file : in file_type; p : in Poly );
  procedure put_line ( file : in file_type; p : in Poly;
                       s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Starts another line after the output of every term.

  -- ON ENTRY :
  --   file     optional file, must be opened for output;
  --   p        Laurent polynomial in several variables;
  --   s        optional array of symbols for the variables.

  procedure Display_Format;

  -- DESCRIPTION :
  --   Displays formatting information.

end Multprec_Complex_Laurentials_io;
