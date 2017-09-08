with text_io;                            use text_io;
with Ada.Calendar;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Black_Box_Solvers is

-- DESCRIPTION :
--   This package contains the solve operations for the blackbox solver.
--   The solvers arise in several flavors:
--   1) silent (true or false) or with output to file (2),
--   2) for ordinary polynomial or Laurent systems (2),
--   3) without multitasking or with multitasking (2),
--   4) in double, double double, or quad double precision (3).
--   The combinations result in 24 Solve procedures.

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List );
  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   without output to file, available for three levels of precision:
  --   standard double, double double, and quad double.

  -- ON INPUT :
  --   p            a polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   sols         solutions found at the end of the paths.

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List );
  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p.

  -- ON INPUT :
  --   file         must be opened for output;
  --   p            a polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   sols         solutions found at the end of the paths.

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List );
  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the Laurent system p,
  --   without multitasking and without intermediate output to file,
  --   in standard double, double double, or quad double precision.

  -- ON INPUT :
  --   p            a Laurent polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   sols         solutions found at the end of the paths.

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List );
  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the Laurent system p,
  --   without multitasking and with output to file,
  --   in standard double, double double, or quad double precision.

  -- ON INPUT :
  --   file         must be opened for output;
  --   p            a Laurent polynomial system.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   sols         solutions found at the end of the paths.

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   using nt tasks, without output to file, and with
  --   standard double, double double, or quad double arithmetic.

  -- ON INPUT :
  --    nt          number of tasks for multithreading, 0 if no multitasking;
  --    p           a polynomial system.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   sols         solutions found at the end of the paths.

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   using nt tasks, with output to file, and with
  --   standard double, double double, or quad double arithmetic.

  -- ON INPUT :
  --   file         must be opened for output;
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a polynomial system.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   sols         solutions found at the end of the paths.

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   with nt tasks and no intermediate output to file,
  --   in standard double, double double, or quad double precision.

  -- ON INPUT :
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a Laurent polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   sols         solutions found at the end of the paths.

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the Laurent polynomial system p,
  --   using nt tasks in standard double, double double, or quad double,
  --   with intermediate output written to file.

  -- ON INPUT :
  --   file         must be opened for output;
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a Laurent polynomial system.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   sols         solutions found at the end of the paths.

end Black_Box_Solvers;
