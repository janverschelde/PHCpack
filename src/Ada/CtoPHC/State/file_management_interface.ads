with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package File_Management_Interface is

-- DESCRIPTION :
--   Provides function to manage input and output files.

  function File_Management_Prompt_Input_File
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for the name of an input file and
  --   opens the file for input.
  --   The verbose level is given by the value of vrblvl.

  function File_Management_Prompt_Output_File
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for the name of an output file and
  --   and makes a file for output.
  --   The verbose level is given by the value of vrblvl.

  function File_Management_Set_Output
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines the output file with the given file name.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       contains the name of the output file;
  --   vrblvl  is the verbose level.

  function File_Management_Close_Output
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Closes the defined output file.
  --   The verbose level is given in vrblvl.

  function File_Management_Write_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes a string to the defined output file.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       contains the string which will be written to file;
  --   vrblvl  is the verbose level.

  function File_Management_Write_Integers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes a sequence of integers to the defined output file.

  -- ON ENTRY :
  --   a       in a[0] is the number of integers in b;
  --   b       the integers to be written to file;
  --   vrblvl  is the verbose level.

  function File_Management_Write_Doubles
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes a sequence of doubles to the defined output file.

  -- ON ENTRY :
  --   a       in a[0] is the number of integers in b;
  --   c       the doubles to be written to file;
  --   vrblvl  is the verbose level.

end File_Management_Interface;
