with text_io;                            use text_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Prompt_for_Solutions is

-- DESCRIPTION :
--   By default, the solutions are on the same file as the system.
--   The procedures in this package scan the file for the solutions
--   and prompt the user in case of failure.

  procedure Scan_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out Standard_Complex_Solutions.Solution_List;
                found : out boolean );
  procedure Scan_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                found : out boolean );
  procedure Scan_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                found : out boolean );

  -- DESCRIPTION :
  --   Attempts to read the solutions from the file, if onfile.

  -- ON ENTRY :
  --   file     opened for input if onfile;
  --   onfile   indicates whether the input system is on file.

  -- ON RETURN :
  --   file     closed if onfile;
  --   sols     solutions if onfile and if found;
  --   found    indicates whether solutions are on file.

  procedure Read_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Read_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Read_Solutions
              ( file : in out file_type; onfile : in boolean;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Attempts to read the solutions from the file, if onfile.
  --   If not onfile, or if the reading from file failed,
  --   then the user is prompted to provide a file for the solutions.

  -- ON ENTRY :
  --   file     opened for input if onfile;
  --   onfile   indicates whether the input system is on file.

  -- ON RETURN :
  --   file     closed for input if onfile;
  --   sols     solutions read from file.

end Prompt_for_Solutions;
