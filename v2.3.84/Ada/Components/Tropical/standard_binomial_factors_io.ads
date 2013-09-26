with text_io;                            use text_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Binomial_Factors;          use Standard_Binomial_Factors;

package Standard_Binomial_Factors_io is

-- DESCRIPTION :
--   Provides output of initial terms of a Puiseux series expansion.

  procedure Initialize_Symbol_Table ( s : in string );
  procedure Initialize_Symbol_Table ( s1,s2 : in string );

  -- DESCRIPTION :
  --   Initializes the symbol table with the symbols in s, s1, s2.

  procedure put ( t : in Term );
  procedure put ( file : in file_type; t : in Term );

  -- DESCRIPTION :
  --   Writes the information in the term as a leading term
  --   of a Puiseux series expansion.

  procedure put ( t : in List_of_Terms );
  procedure put ( file : in file_type; t : in List_of_Terms );

  -- DESCRIPTION :
  --   Writes all terms in the list to standard output or to file.

end Standard_Binomial_Factors_io;
