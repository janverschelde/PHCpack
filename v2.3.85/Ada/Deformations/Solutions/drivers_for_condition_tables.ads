with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Drivers_for_Condition_Tables is

-- DESCRIPTION :
--   The procedures defined by this package are some basic utilities
--   and main interactive drivers for a condition report on the list
--   of computed isolated solutions of a polynomial system.

  procedure Read_and_Compute_Condition_Tables;

  -- DESCRIPTION :
  --   This is the first "Main" program to compute a condition table.
  --   It may crash for huge lists due to insufficient memory.

  procedure Read_Dimensions
               ( file : in file_type; len,dim : out natural32;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Reads the length of a solution list and the dimension
  --   of its solution vectors from file and reports failure.

  procedure Prompt_to_Scan_Banner
              ( infile : in file_type; bannered : out boolean;
                fail : out boolean );

  -- DESCRIPTION :
  --   Prompts the user whether the solutions on the infile are preceeded
  --   by a system.  The fail is true if the banner was not found.

  -- ON RETURN :
  --   bannered is true if the user indicated that the solutions
  --            are preceeded by a banner, false otherwise;
  --   fail     true if bannered and if the banner was not found;
  --            false otherwise.

  procedure Scan_Banner_Dimensions
              ( infile : in file_type; len,dim : out natural32;
                bannered,fail : out boolean );

  -- DESCRIPTION :
  --   Prompts the user, asking whether the solutions on file are
  --   preceeded by a banner, before reading the solution dimensions.

  -- ON ENTRY :
  --   infile   input file, must be opened for reading.

  -- ON RETURN :
  --   len      length of the solution list on file, if not fail;
  --   dim      dimensions of the solutions on file, if not fail;
  --   bannered is true if the user indicated that the solutions are
  --            proceeded by a banner, false otherwise;
  --   fail     true if the length and dimensions could not be read.

  procedure Interactive_Driver_to_Scan_Solution_Lists;

  -- DESCRIPTION :
  --   This is the main interactive driver to scan a list of solutions
  --   and compute a condition table.  The user is prompted to provide
  --   input and output file.  The important key feature is that the
  --   list is not stored entirely at once into main internal memory.

  procedure Main_Driver_to_Scan_Solution_Lists
              ( infilename,outfilename : in string );
 
  -- DESCRIPTION :
  --   This is the driver as called by the main validation module in phc.
  --   The strings given on entry are respectively the names of the input
  --   and output files.  If the string for the output file is empty,
  --   then output will be written to screen.  If the string for the 
  --   input file is empty, then the user will be prompted to provide
  --   a name for the input file.

end Drivers_for_Condition_Tables;
