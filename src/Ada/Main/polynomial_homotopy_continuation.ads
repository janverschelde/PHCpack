with text_io;                            use text_io;
with Ada.Calendar;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with String_Splitters;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Polynomial_Homotopy_Continuation is

-- DESCRIPTION :
--   Defines phc when call without any options, or as phc -t.

  procedure Display_all_Options ( nt : in natural32 );

  -- DESCRIPTION :
  --   Displays an overview of all options on screen.

  procedure Main_Polynomial_Solver 
              ( file : in file_type; nt : in natural32; p : in Poly_Sys;
                ls : in String_Splitters.Link_to_Array_of_Strings;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs interactively through the stages to solve p,
  --   writing output to file.  The array of strings contains the
  --   original string representation of p, and may be needed if
  --   the working precision is increased.
  --   The value for the verbose level is given by vrb.

  procedure Main_Laurent_Solver 
              ( file : in file_type; nt : in natural32; p : in Laur_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   This is the main interactive solver for Laurent systems,
  --   primarily guiding through polyhedral homotopies.

  procedure Start_Main
              ( start_moment : in Ada.Calendar.Time; nt : in natural32;
                outfilename : in string; q : in Laur_Sys;
                ls : in String_Splitters.Link_to_Array_of_Strings;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Starts the running through the stages,
  --   converting q to an ordinary polynomial system
  --   if no negative exponents are present.

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs phc in 'full mode' without specific options,
  --   other than the number of tasks.

  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used to track the paths;
  --   infilename     the name of the input file;
  --   outfilename    the name of the output file;
  --   verbose        the verbose level.

end Polynomial_Homotopy_Continuation;
