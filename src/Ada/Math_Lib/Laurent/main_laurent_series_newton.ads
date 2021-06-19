with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Main_Laurent_Series_Newton is

-- DESCRIPTION :
--   The main procedures for phc -u to run Newton's method on Laurent series.
--   Only double precision is currently supported.

  procedure Start_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Starts Newton's method at constant terms.
  --   Reads the polynomial system with start solutions from the file
  --   defined by infilename, or prompts for a system and solutions.
  --   The name of the output file is in outfilename.

  procedure Start_at_Series
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Starts Newton's method at some series.
  --   Reads the polynomial system from the file defined by infilename,
  --   or prompts for a system.

  procedure Run_Laurent_Series_Newton
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method at Laurent series,
  --   after prompting whether to start at a constant or at a series.
  --   The file name for the input system may be given in infilename.
  --   The file name for the output may be given in outfilename.
  --   If empty strings are passed, then the user is prompted.
  --   The value of the verbose level is in vrb.

end Main_Laurent_Series_Newton;
