with Ada.Calendar;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package Black_Box_Square_Solvers is

-- DESCRIPTION :
--   A square system has as many equations as unknowns
--   and are the input to the solve procedures below.

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                    deflate,append_sols : in boolean;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   A polynomial system with as many equations as unknowns is square.
  --   This procedure solves a square polynomial system p,
  --   using double, double double, or quad double arithmetic.
  --
  -- REQUIRED :
  --   Must be called after the special cases (one single equation
  --   and a linear system) have been dealt with.

  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used to track the paths;
  --   infilename     the name of the input file;
  --   outfilename    the name of the output file;
  --   start_moment   clock time when phc was started;
  --   p              polynomial system to be solved.
  --   deflate        if not deflate, then not deflation will be applied;
  --   append_sols    true if solutions need to be appended to input file;
  --   verbose        the verbose level.

  -- ON RETURN :
  --   p              system may be scaled or reduced.

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   A polynomial system with as many equations as unknowns is square.
  --   This procedure solves a square Laurent polynomial system p,
  --   with double, double double, or quad double arithmetic.
  --
  -- REQUIRED :
  --   Must be called after the special cases (one single equation
  --   and a linear system) have been dealt with.

  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used to track the paths;
  --   infilename     the name of the input file;
  --   outfilename    the name of the output file;
  --   start_moment   clock time when phc was started;
  --   p              polynomial system to be solved.
  --   append_sols    true if solutions need to be appended to input file;
  --   verbose        the verbose level.

  -- ON RETURN :
  --   p              system may be scaled or reduced.

end Black_Box_Square_Solvers;
