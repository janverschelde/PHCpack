with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with DoblDobl_Complex_Solutions;        use DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Continuation_Data;        use DoblDobl_Continuation_Data;

package Drivers_to_Track_DoblDobl_Paths is

-- DESCRIPTION :
--   This package exports path trackers in double double arithmetic
--   with incremental read and write of the solution lists.

  generic

    with procedure Get_Next ( n : in natural32; ls : out Link_to_Solution;
                              ind : out natural32; continue : out boolean );
    -- DESCRIPTION :
    --   Returns in ls a new start solution of dimension n with index ind.
    --   If the ls on return is null, then the tracking stops.
    --   If continue is false, then Get_Next will not be called again.

  procedure Track ( file : in file_type; report : in boolean;
                    p,q : in Poly_Sys; len,dim : in natural32;
                    f : access procedure ( s : in Solu_Info ) := null );

  -- DESCRIPTION :
  --   Tracks paths using a homotopy between q and p.  Whenever a new start
  --   solution is needed, the procedure Get_Next is called.

  -- ON ENTRY :
  --   file     to write the end points of the solution paths;
  --   report   if true, then user can monitor progress of solver on screen,
  --            as one line will be written for every start solution;
  --   p        target system;
  --   q        start system;
  --   len      total number of paths to be tracked;
  --   dim      dimension of the solution vectors;
  --   f        the callback function during path tracking.

  procedure Cheater_Track
              ( outfile,solfile : in file_type; report : in boolean;
                p,q : in Poly_Sys; len,dim : in natural32;
                f : access procedure ( s : in Solu_Info ) := null );

  -- DESCRIPTION :
  --   Tracks paths, starting at solutions of q to solve p,
  --   reading and writing solutions one at a time.

  -- REQUIRED :
  --   The solution file solfile is positioned appropriately at
  --   the next start solution of q.

  -- ON ENTRY :
  --   outfile  to write the end points of the solution paths;
  --   solfile  is the file where the list of solutions are;
  --   report   if true, then user can monitor progress of solver on screen,
  --            as one line will be written for every start solution;
  --   p        target system;
  --   q        start system;
  --   len      total number of paths to be tracked;
  --   dim      dimension of the solution vectors;
  --   f        the callback function during path tracking.

  procedure Total_Degree_Track
              ( file : in file_type; report : in boolean;
                p,q : in Poly_Sys; len,dim : in natural32;
                f : access procedure ( s : in Solu_Info ) := null );

  -- DESCRIPTION :
  --   Tracks paths, starting at solutions of q to solve p,
  --   computing a new root of q, tracking the path and writing
  --   the result to file, one solution at a time.

  -- REQUIRED :
  --   The start system q is based on the total degree.

  -- ON ENTRY :
  --   file     to write the end points of the solution paths;
  --   report   if true, then user can monitor progress of solver on screen,
  --            as one line will be written for every start solution;
  --   p        target system;
  --   q        start system;
  --   len      total number of paths to be tracked;
  --   dim      dimension of the solution vectors;
  --   f        the callback function during path tracking.

  procedure Get_Lex_Linear_Product_Start_Solution
              ( report : in boolean; n : in natural32;
                d,cp : in Standard_Natural_Vectors.Vector;
                ind : in natural32; ndp : in natural32;
                tol : in double_float;
                ls : out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Solves one linear system defined by the current index ind,
  --   for one solution of a linear-product start system.

  -- ON ENTRY :
  --   report   flag to indicate if monitoring on screen is wanted;
  --   n        number of equations and variables;
  --   d        degrees of the equations in the start system;
  --   cp       consecutive products of the degrees in d;
  --   ind      index to the current solution;
  --   ndp      number of decimal places needed in cnt,
  --            the value of ndp matters only if report is true;
  --   tol      tolerance on the condition number of a good solution.

  -- ON RETURN :
  --   ls       the next solution of the linear-product start system;
  --   fail     true if there is no next start solution.

  procedure Next_Lex_Linear_Product_Start_Solution
              ( report : in boolean; n : in natural32;
                d,cp : in Standard_Natural_Vectors.Vector;
                cnt : in out natural32; len,ndp : in natural32;
                tol : in double_float;
                ls : out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Computes the next solution of a linear-product start system,
  --   executing a loop starting at index cnt+1,
  --   using a plain lexicographic enumeration.

  -- ON ENTRY :
  --   report   flag to indicate if monitoring on screen is wanted;
  --   n        number of equations and variables;
  --   d        degrees of the equations in the start system;
  --   cp       consecutive products of the degrees in d;
  --   cnt      number of solutions already computed;
  --   len      total degree of the start system;
  --   ndp      number of decimal places needed in cnt,
  --            the value of ndp matters only if report is true;
  --   tol      tolerance on the condition number of a good solution.

  -- ON RETURN :
  --   cnt      incremented counter on number of linear systems solved;
  --   ls       the next solution of the linear-product start system;
  --   fail     true if there is no next start solution.

  procedure Get_Next_Linear_Product_Start_Solution
              ( report : in boolean; n : in natural32;
                s : in out Standard_Natural_Vectors.Vector;
                d : in Standard_Natural_Vectors.Vector;
                cnt : in out natural32; len,ndp : in natural32;
                tol : in double_float;
                ls : out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Computes the next solution of a linear-product start system.

  -- ON ENTRY :
  --   report   flag to indicate if monitoring on screen is wanted;
  --   n        number of equations and variables;
  --   s        position vector for the current solution;
  --   d        degrees of the equations in the start system;
  --   cnt      number of solutions already computed,
  --            must be zero for the first solution;
  --   len      total degree of the start system;
  --   ndp      number of decimal places needed in cnt,
  --            the value of ndp matters only if report is true;
  --   tol      tolerance on the condition number of a good solution.

  -- ON RETURN :
  --   s        new position vector for a solution if not fail;
  --   cnt      incremented counter on number of linear systems solved;
  --   ls       the next solution of the linear-product start system;
  --   fail     true if there is no next start solution.

  procedure Get_Next_Linear_Product_Start_Solution
              ( report : in boolean; n : in natural32;
                cnt : in out natural32; len,ndp : in natural32;
                tol : in double_float;
                ls : out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Computes the next solution of a linear-product start system.
  --   This version of Get_Next has memory of the "s" and "d".

  -- ON ENTRY :
  --   report   flag to indicate if monitoring on screen is wanted;
  --   n        number of equations and variables;
  --   cnt      number of solutions already computed,
  --            must be zero for the first solution;
  --   len      total degree of the start system;
  --   ndp      number of decimal places needed in cnt,
  --            the value of ndp matters only if report is true;
  --   tol      tolerance on the condition number of a good solution.

  -- ON RETURN :
  --   cnt      incremented counter on number of linear systems solved;
  --   ls       the next solution of the linear-product start system;
  --   fail     true if there is no next start solution.

  procedure Linear_Product_Track
              ( file : in file_type; report : in boolean;
                p,q : in Poly_Sys; len,dim : in natural32;
                f : access procedure ( s : in Solu_Info ) := null );

  -- DESCRIPTION :
  --   Tracks paths, starting at solutions of q to solve p,
  --   computing a new root of q, tracking the path and writing
  --   the result to file, one solution at a time.

  -- REQUIRED :
  --   The start system q is a linear-product start system,
  --   properly stored in standard_linear_product_system.

  -- ON ENTRY :
  --   file     to write the end points of the solution paths;
  --   report   if true, then user can monitor progress of solver on screen,
  --            as one line will be written for every start solution;
  --   p        target system;
  --   q        start system;
  --   len      total number of paths to be tracked,
  --            this may be initialized with the total degree;
  --   dim      dimension of the solution vectors;
  --   f        the callback function during path tracking.

  procedure Read_Target_System
              ( file : in file_type; p : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Reads the target system from the file, which must be opened for input.
  --   Exceptions are reported to the user and raised.
  --   On return in p is the target system.

  procedure Read_Start_System
              ( file : in file_type; p : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Reads the start system from the file, which much be opened for input.
  --   Exceptions are reported to the user and raised.
  --   On return in p is the start system.

  procedure Read_Systems_and_Track 
              ( pfile,qfile,ofile : in file_type;
                f : access procedure ( s : in Solu_Info ) := null );

  -- DESCRIPTION :
  --   Reads the target and start systems from file and calls the
  --   other Read_Systems_and_Track below.

  -- REQUIRED :
  --   pfile and qfile are already opened for input;
  --   the file ofile has been newly created for output.

  -- ON ENTRY :
  --   pfile    file with the target system;
  --   qfile    file with the start system (and start solutions);
  --   ofile    output file for diagnostics and results;
  --   f        the callback function during path tracking.

  procedure Read_Systems_and_Track
              ( target : in Poly_Sys; qfile,ofile : in file_type;
                kind : in character;
                f : access procedure ( s : in Solu_Info ) := null );

  -- DESCRIPTION :
  --   Depending on the type of start system, the appropriate version
  --   of Track is launched.

  -- REQUIRED :
  --   qfile is open for input and ofile is ready for output.

  -- ON ENTRY :
  --   target   polynomial system that is to be solved;
  --   qfile    file with the start system (and start solutions);
  --   ofile    output file for diagnostics and results;
  --   kind     kind of start system that will be used, must be
  --            '1' (total degree) or '3' (cheater);
  --   f        the callback function during path tracking.

end Drivers_to_Track_DoblDobl_Paths;
