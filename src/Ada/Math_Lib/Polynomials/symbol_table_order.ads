with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Symbol_Table_Order is

-- DESCRIPTION :
--   Defines phc -o.

  procedure Main ( infilename,outfilename : in string;
                   vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads a polynomial system from the input file
  --   and writes the symbols in the symbol table to the output file.
  --   The value of the verbose level is in vrblvl.
  --   When the input string are empty, the user is prompted for
  --   the input and output files.

end Symbol_Table_Order;
