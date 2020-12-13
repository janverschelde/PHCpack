with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;

package Black_Box_Polyhedral_Solvers is

-- DESCRIPTION :
--   The solvers in this package assume square polynomial systems
--   and applied polyhedral homotopies to solve those systems.

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

  -- DESCRIPTION :
  --   Applies polyhedral homotopies in blackbox mode to solve p,
  --   without output to file.  Either the output of the root counter 
  --   is written to standard output or to a string.

  -- ON INPUT :
  --   p            a polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen;
  --   deflate      if not deflate, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   rocos        string with the root count information,
  --                displayed in the format as when silent is false
  --                in the first Polyhedral_Solve;
  --   gamma        complex gamma constant used in homotopy;
  --   q            start system;
  --   qsols        start solutions;
  --   sols         solutions found at the end of the paths.

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 );
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

  -- DESCRIPTION :
  --   Applies polyhedral homotopies in blackbox mode to solve p,
  --   without output to file.  Either the output of the root counter 
  --   is written to standard output or to a string.

  -- ON INPUT :
  --   nt           the number of tasks,
  --                if nt < 2, then Solve without nt is called;
  --   p            a polynomial system;
  --   silent       if true, then the computed root counts will not be shown,
  --                if false, then the user will see the computed root counts
  --                displayed on screen;
  --   deflate      if not deflate, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   rocos        string with the root count information,
  --                displayed in the format as when silent is false
  --                in the first Polyhedral_Solve;
  --   gamma        random complex gamma constant used;
  --   q            start system;
  --   qsols        start solutions;
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

  -- DESCRIPTION :
  --   Applies polyhedral homotopies to solve p with output to file.

  -- ON INPUT :
  --   file         must be opened for output;
  --   nt           the number of tasks,
  --                if nt < 2, then Solve without nt is called;
  --   p            a polynomial system;
  --   deflate      if false, then no deflation will be applied;
  --   verbose      the verbose level.

  -- ON RETURN :
  --   rc           root count used in the homotopy to solve p;
  --   gamma        random complex gamma constant used;
  --   q            start system;
  --   qsols        start solutions;
  --   sols         solutions found at the end of the paths.

end Black_Box_Polyhedral_Solvers;
