with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Laur_Systems;

package Main_Laurent_Series_Newton is

-- DESCRIPTION :
--   The main procedures for phc -u to run Newton's method on Laurent series.
--   Only double precision is currently supported.

  procedure Run_Regular_Newton
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys; 
                x : in Standard_Complex_Vectors.Vector;
                tdx,deg : in integer32; verbose : in boolean );

  -- DESCRIPTION :
  --   Runs Newton's method on the system p,
  --   with series starting at the constants in x.

  -- ON ENTRY :
  --   file     must be open for output;
  --   p        system with one parameter t;
  --   x        solution for t = 0;
  --   tdx      index of t as a variable in p;
  --   deg      truncation degree of the series;
  --   verbose  is the verbose flag.

  procedure Run_Singular_Newton
              ( file : in file_type;
                p,x : in Standard_Complex_Laur_Systems.Laur_Sys; 
                tdx,deg : in integer32; verbose : in boolean );

  -- DESCRIPTION :
  --   Runs Newton's method on the system p,
  --   with series starting at the constants in x.

  -- ON ENTRY :
  --   file     must be open for output;
  --   p        system with one parameter t;
  --   x        initial terms of a power series;
  --   tdx      index of t as a variable in p;
  --   deg      truncation degree of the series;
  --   verbose  is the verbose flag.

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
