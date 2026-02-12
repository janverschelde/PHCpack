with Ada.Text_IO;                        use Ada.Text_IO;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;

package Standard_Complex_Laur_Readers is

-- DESCRIPTION :
--   The procedures in this packages are the auxiliaries for the
--   get procedures in the standard_complex_laurentials_io package.

  procedure Read_Term ( file : in file_type; bc : in out integer32;
                        char : in out character;
                        n : in natural32; termp : in out Poly );
  -- DESCRIPTION :
  --   Reads a term from file, char is the first character of the term.

  -- ON ENTRY :
  --   file     must be opened for input;
  --   bc       counts #open - #closed brackets;
  --   char     first character of the term;
  --   n        number of variables;
  --   termp    accumulating polynomial for the term.

  -- ON RETURN :
  --   bc       updated bracket counter;
  --   char     first character following the term;
  --   termp    resulting term read from file.

  procedure Read_Factor ( file : in file_type; bc : in out integer32;
                          char : in out character; n : in natural32;
                          d : in out Degrees; pb : in out Poly );
  -- DESCRIPTION :
  --   Reads a factor from file, char is the first character of the factor.

  -- DESCRIPTION :
  --   Reads a factor from file, char is the first character of the factor.

  -- ON ENTRY :
  --   file     must be opened for input;
  --   bc       counts #open - #closed brackets;
  --   char     first character of the factor;
  --   n        number of variables;
  --   d        accumulates the degrees of the factor;
  --   pb       accumulating polynomial for the factor.

  -- ON RETURN :
  --   bc       updated bracket counter;
  --   char     first character following the factor;
  --   pb       resulting factor read from file.

  procedure Read_Power_Factor
              ( file : in file_type; char : in out character;
                p : in out Poly );

  -- DESCRIPTION :
  --   Reads the exponent of a factor stored in p and returns
  --   in p the polynomial p raised to the read exponent.

  -- ON ENTRY :
  --   file     must be opened for input;
  --   char     equals the exponentiation symbol '^';
  --   p        current factor.

  -- ON RETURN :
  --   char     character following the exponent;
  --   p        p raised to the read exponent.
 
  procedure Read_Polynomial
              ( file : in file_type; bc : in out integer32;
                p : out Poly );

  -- DESCRIPTION :
  --   Reads a polynomial from file, returned in p, raising BAD_BRACKET
  --   if the bracket counter bc is nonzero when the delimiter is encountered.

end Standard_Complex_Laur_Readers;
