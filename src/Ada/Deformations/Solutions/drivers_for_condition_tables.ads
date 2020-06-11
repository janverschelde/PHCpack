with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Drivers_for_Condition_Tables is

-- DESCRIPTION :
--   The procedures defined by this package are some basic utilities
--   and main interactive drivers for a condition report on the list
--   of computed isolated solutions of a polynomial system.

  procedure Standard_Read_and_Compute_Condition_Tables;
  procedure DoblDobl_Read_and_Compute_Condition_Tables;
  procedure QuadDobl_Read_and_Compute_Condition_Tables;

  -- DESCRIPTION :
  --   This is the first "Main" program to compute a condition table,
  --   for a list of solutions in standard double precision (Standard),
  --   double double (DoblDobl), or quad double (QuadDobl) precision.
  --   Because this procedure reads the entire solution list into main
  --   memory, it may crash for huge lists due to insufficient memory.

  procedure Interactive_Driver_to_Scan_Solution_Lists;

  -- DESCRIPTION :
  --   This is the main interactive driver to scan a list of solutions
  --   and compute a condition table.  The user is prompted to provide
  --   input and output file.  The important key feature is that the
  --   list is not stored entirely at once into main internal memory.

  procedure Main_Driver_to_Scan_Solution_Lists
              ( infilename,outfilename : in string;
                verbose : in integer32 := 0 );
 
  -- DESCRIPTION :
  --   This is the driver as called by the main validation module in phc.
  --   The strings given on entry are respectively the names of the input
  --   and output files.  If the string for the output file is empty,
  --   then output will be written to screen.  If the string for the 
  --   input file is empty, then the user will be prompted to provide
  --   a name for the input file.  The last parameter is the verbose level
  --   and if > 0, then one line is written to screen.

end Drivers_for_Condition_Tables;
