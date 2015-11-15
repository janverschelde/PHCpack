with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Symbol_Table;                       use Symbol_Table;

package Write_Factors is

-- DESCRIPTION :
--   The procedures in this package define the output of one single factor,
--   that is of the form x(i)^d for one single variable x(i) to the power d.
--   These procedure are shared for all output of all polynomials,
--   regardless the type of the coefficients.

  procedure Write_Factor ( file : in file_type; d : in natural32;
                           sb : in Symbol; pow : in power );

  -- DESCRIPTION :
  --   Writes the factor corresponding with the i-th unknown on file,
  --   using the symbol in sb to represent the i-th variable.

  procedure Write_Factor ( file : in file_type; d,i : in natural32;
                           standard : in boolean; pow : in power );

  -- DESCRIPTION :
  --   Writes the factor corresponding with the i-th variable on file.
  --   If standard, then 'x' will be used as the symbol,
  --   otherwise the symbol for the i-th variable will be retrieved
  --   from the symbol table.

end Write_Factors;
