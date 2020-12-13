with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Black_Box_Solvers is

-- DESCRIPTION :
--   This package contains the solve operations for the blackbox solver.
--   The solvers arise in several flavors:
--   1) silent (true or false) or with output to file,
--   2) for ordinary polynomial or Laurent systems,
--   3) without multitasking or with multitasking,
--   4) in double, double double, or quad double precision,
--   5) with root counts to string or not,
--   6) with the return or not of start system and start solutions.
--   The combinations result in 72 Solve procedures.

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   without output to file, available for three levels of precision:
  --   standard double, double double, and quad double.

  -- ON INPUT :
  --   p            a polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen;
  --   deflate      if not deflate, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   without output to file, available for three levels of precision:
  --   standard double, double double, and quad double.

  -- ON INPUT :
  --   p            a polynomial system;
  --   deflate      if false, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   rocos        string with the root count information,
  --                displayed in the format as when silent is false
  --                in the other above solve procedures,
  --                rocos is null if p is one of the special cases!;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p.

  -- ON INPUT :
  --   file         must be opened for output;
  --   p            a polynomial system;
  --   deflate      if false, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the Laurent system p,
  --   without multitasking and without intermediate output to file,
  --   in standard double, double double, or quad double precision.

  -- ON INPUT :
  --   p            a Laurent polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the Laurent system p,
  --   without multitasking and without intermediate output to file,
  --   in standard double, double double, or quad double precision.

  -- ON INPUT :
  --   p            a Laurent polynomial system;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   rocos        information about the root counts, in the same format
  --                as the above solve procedures with false for silent,
  --                rocos is null if p is one of the special cases!;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the Laurent system p,
  --   without multitasking and with output to file,
  --   in standard double, double double, or quad double precision.

  -- ON INPUT :
  --   file         must be opened for output;
  --   p            a Laurent polynomial system;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean; rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   using nt tasks, without output to file, and with
  --   standard double, double double, or quad double arithmetic.

  -- ON INPUT :
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a polynomial system.
  --   silent       if not silent, then root counting information will be
  --                written the standard output;
  --   deflate      if false, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   using nt tasks, without output to file, and with
  --   standard double, double double, or quad double arithmetic.

  -- ON INPUT :
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a polynomial system;
  --   deflate      if false, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   rocos        string with the root count information,
  --                displayed in the format as when silent is false
  --                in the other above solve procedures,
  --                rocos is null if p is one of the special cases!;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   using nt tasks, with output to file, and with
  --   standard double, double double, or quad double arithmetic.

  -- ON INPUT :
  --   file         must be opened for output;
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a polynomial system;
  --   deflate      if false, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system (if the system was not special);
  --   qsols        start solutions (if the system was not special);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   with nt tasks and no intermediate output to file,
  --   in standard double, double double, or quad double precision.

  -- ON INPUT :
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a Laurent polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system (if p is not a special case system);
  --   qsols        start solutions (if p is not a special case system);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the polynomial system p,
  --   with nt tasks and no intermediate output to file,
  --   in standard double, double double, or quad double precision.

  -- ON INPUT :
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a Laurent polynomial system;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   rocos        string with the root count information,
  --                displayed in the format as when silent is false
  --                in the other above solve procedures,
  --                rocos is null if p is one of the special cases!;
  --   gamma        random complex gamma constant used;
  --   q            start system (if p is not a special case system);
  --   qsols        start solutions (if p is not a special case system);
  --   sols         solutions found at the end of the paths.

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the blackbox solver to solve the Laurent polynomial system p,
  --   using nt tasks in standard double, double double, or quad double,
  --   with intermediate output written to file.

  -- ON INPUT :
  --   file         must be opened for output;
  --   nt           number of tasks for multithreading, 0 if no multitasking;
  --   p            a Laurent polynomial system;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system (if p is not a special case system);
  --   qsols        start solutions (if p is not a special case system);
  --   sols         solutions found at the end of the paths.

end Black_Box_Solvers;
