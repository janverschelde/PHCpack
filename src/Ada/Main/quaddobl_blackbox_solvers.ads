with Ada.Calendar;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package QuadDobl_BlackBox_Solvers is

-- DESCRIPTION :
--   Defines what phc -b4 does.
--   For Laurent binomial systems (the genuine ones with negative powers),
--   a stable mixed volume or an affine solution set does not make sense.

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    append_sols : in boolean; v : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the blackbox solver for a polynomial system.

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                    append_sols : in boolean; v : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the blackbox solver for a Laurent polynomial system.

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This is the main interactive driver for the homotopy continuation
  --   package for the blackbox solution of polynomial systems,
  --   with computations done with quad double arithmetic.
  --   This routine is executed with the option -b4 of phc.

  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used to track the paths;
  --   infilename     the name of the input file;
  --   outfilename    the name of the output file;
  --   verbose        the verbose level.

end QuadDobl_BlackBox_Solvers;
