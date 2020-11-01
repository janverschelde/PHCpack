with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Main_Dictionary_Solutions is

-- DESCRIPTION :
--   Defines the format conversions of solution lists from and into
--   Python dictionaries.  Main is executed by phc -x.

  procedure Scan_Solutions
               ( filename : in string; sols : in out Solution_List;
                 solsonfile : out boolean; python_format : out boolean );

  -- DESCRIPTION :
  --   Checks whether the given file name corresponds to a file with
  --   the solutions in a correct format.  If so, then solsonfile is
  --   true on return and sols contain the solutions.
  --   When the solutions are in Python format, then python_format is true.

  procedure Write ( file : in file_type; python_format : in boolean;
                    sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the solutions to the file.

  -- REQUIRED : file is opened in the right mode.

  procedure Write_Solutions
              ( outfilename : in string; python_format : in boolean;
                sols : in Solution_List );

  -- DESCRIPTION :
  --   If the file with name outfilename does not exist,
  --   then it is created and used to write the solutions on;
  --   otherwise, the solutions are appended to the file with
  --   name outfilename.

  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This routine converts output solutions into a Python dictionary
  --   as called by phc -x.
  --   The first two arguments are the names of the input and output file.
  --   The last argument is the verbose level.

  -- FUNCTIONALITY :
  --   The program strips the output file of phc (with name infilename)
  --   of all intermediate output, scans that file till "THE SOLUTIONS"
  --   is reached and then converts the solutions into a Python dictionary.
  --   If the outfilename is the empty string, then the output will
  --   be written to the screen, otherwise to the file name.

end Main_Dictionary_Solutions;
