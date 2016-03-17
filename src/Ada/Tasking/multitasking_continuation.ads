with text_io;                           use text_io;
with String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Multitasking_Continuation is

-- DESCRIPTION :
--   Offers routines to track all paths using multitasking.

  procedure Silent_Path_Tracker
               ( ls : in Standard_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );
  procedure Silent_Path_Tracker
               ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );
  procedure Silent_Path_Tracker
               ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Wrapper to track one path in standard double, double double,
  --   or quad double precision, starting at the solution in ls.
  --   There is no intermediate output of any kind.
  --   If nbq differs from the zero default and equals the number of
  --   equations, then Gauss-Newton correctors are applied.

  procedure Silent_Path_Tracker
               ( id,nb : in integer32;
                 ls : in Standard_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );
  procedure Silent_Path_Tracker
               ( id,nb : in integer32;
                 ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );
  procedure Silent_Path_Tracker
               ( id,nb : in integer32;
                 ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Task with identification number id reports the receipt of
  --   solution with number nb, with data in ls.
  --   With id and nb given,  one line is written to screen allowing the 
  --   user to monitor the progress of the computations.
  --   If nbq differs from the zero default and equals the number of
  --   equations, then Gauss-Newton correctors are applied.

  procedure Silent_Laurent_Path_Tracker
               ( ls : in Standard_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );
  procedure Silent_Laurent_Path_Tracker
               ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );
  procedure Silent_Laurent_Path_Tracker
               ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Wrapper to track one path defined by a Laurent homotopy in
  --   standard double, double double, or quad double precision.

  procedure Silent_Laurent_Path_Tracker
               ( id,nb : in integer32;
                 ls : in Standard_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );
  procedure Silent_Laurent_Path_Tracker
               ( id,nb : in integer32;
                 ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );
  procedure Silent_Laurent_Path_Tracker
               ( id,nb : in integer32;
                 ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Path trackers for Laurent polynomial systems.

  procedure Silent_Multitasking_Path_Tracker
               ( sols : in out Standard_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );
  procedure Silent_Multitasking_Path_Tracker
               ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );
  procedure Silent_Multitasking_Path_Tracker
               ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   n tasks will track paths in standard double, double double,
  --   or quad double precision, paths defined by a polynomial homotopy.
  --   If nbq differs from the zero default value and equals the number
  --   of equations in an overdetermined homotopy, then the Gauss-Newton
  --   correctors are applied to track the paths.
  --   No intermediate output is written by this silent version.

  procedure Reporting_Multitasking_Path_Tracker
               ( sols : in out Standard_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );
  procedure Reporting_Multitasking_Path_Tracker
               ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );
  procedure Reporting_Multitasking_Path_Tracker
               ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   n tasks will track paths in standard double, double double,
  --   or quad double precision, paths defined by a polynomial homotopy.
  --   If nbq differs from the zero default value and equals the number
  --   of equations in an overdetermined homotopy, then the Gauss-Newton
  --   correctors are applied to track the paths.
  --   Intermediate output is written to screen by this reporting version:
  --   the user can monitor the progress of the computations.

  procedure Silent_Multitasking_Laurent_Path_Tracker
               ( sols : in out Standard_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );
  procedure Silent_Multitasking_Laurent_Path_Tracker
               ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );
  procedure Silent_Multitasking_Laurent_Path_Tracker
               ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   n tasks will track paths in standard double, double double,
  --   or quad double precision, paths defined by a Laurent homotopy.
  --   If nbq differs from the zero default value and equals the number
  --   of equations in an overdetermined homotopy, then the Gauss-Newton
  --   correctors are applied to track the paths.
  --   No intermediate output is written by this silent version.

  procedure Reporting_Multitasking_Laurent_Path_Tracker
               ( sols : in out Standard_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );
  procedure Reporting_Multitasking_Laurent_Path_Tracker
               ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );
  procedure Reporting_Multitasking_Laurent_Path_Tracker
               ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 n : in integer32; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   n tasks will track the paths starting at sols.
  --   The silent version stays mute, while the reporting version
  --   allows the user to monitor the progress of the computations.

  -- REQUIRED :
  --   {Standard,DoblDobl,QuadDobl}_Homotopy has been initialized properly,
  --   or Standard_Laurent_Homotopy for the *Laurent* trackers.

  -- ON ENTRY :
  --   sols      start solutions.

  -- ON RETURN :
  --   sols      solutions at the end of the paths. 

  procedure Driver_to_Path_Tracker
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 prclvl : in natural32;
                 ls : in String_Splitters.Link_to_Array_of_Strings;
                 n,nbequ,nbvar : in integer32 );
  procedure Driver_to_Path_Tracker
               ( file : in file_type;
                 p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 prclvl : in natural32;
                 ls : in String_Splitters.Link_to_Array_of_Strings;
                 n,nbequ,nbvar : in integer32 );

  -- DESCRIPTION :
  --   Driver to Multitasking_Path_Tracker, for mainpoco.

  -- ON ENTRY :
  --   file      the output file;
  --   p         target polynomial system to be solved;
  --   prclvl    preset level of precision, 1, 2, 4 respectively for 
  --             double, double double, or quad double precision;
  --   ls        string representation of the polynomials in p,
  --             for conversion to double double or quad double precision;
  --   n         number of tasks;
  --   nbequ     the number of equations;
  --   nbvar     the number of variables.

end Multitasking_Continuation;
