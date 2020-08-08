with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Wrapped_Path_Trackers is

-- DESCRIPTION :
--   This package collects wrappers to instantiate a polynomial homotopy
--   and to call the path trackers starting from one solution to track
--   many solution paths in double, double double, and quad double precision.
--   The homotopy continuation parameter is typically the last variable
--   in the polynomial system.

  procedure Set_Parameters ( file : in file_type; report : out boolean );

  -- DESCRIPTION :
  --   Interactive determination of the continuation and output parameters.
  --   The values for the continuation parameters and output code are
  --   written to the file.
  --   The output parameter report is true if the output level was nonzero,
  --   i.e.: if the path trackers are reporting.

-- TRACKING ONE PATH WITHOUT OUTPUT :

  procedure Run ( n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 ); 
  procedure Run ( n : in integer32;
                  h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                  xt : in out DoblDobl_Complex_Vectors.Vector;
                  sol : out DoblDobl_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 ); 
  procedure Run ( n : in integer32;
                  h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                  xt : in out QuadDobl_Complex_Vectors.Vector;
                  sol : out QuadDobl_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 ); 

  -- DESCRIPTION :
  --   Tracks one path starting at the solution in xt,
  --   as defined by the homotopy h, without output,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   n        number of variables in the ambient space;
  --   h        homotopy in n+1 variables;
  --   xt       start solution with its last component equal to zero,
  --            satisfies the homotopy h (upto tolerance);
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   xt       solution at the end of the path, tracked to the
  --            last component of xt to be equal to one;
  --   sol      standard representation of the solution.

-- TRACKING ONE PATH WITH OUTPUT TO FILE :

  procedure Run ( file : in file_type; n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 ); 
  procedure Run ( file : in file_type; n : in integer32;
                  h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                  xt : in out DoblDobl_Complex_Vectors.Vector;
                  sol : out DoblDobl_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 ); 
  procedure Run ( file : in file_type; n : in integer32;
                  h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                  xt : in out QuadDobl_Complex_Vectors.Vector;
                  sol : out QuadDobl_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 ); 

  -- DESCRIPTION :
  --   Tracks one path starting at the solution in xt,
  --   as defined by the homotopy h, with intermediate output,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics,
  --            if omitted, then there is no intermediate output;
  --   n        number of variables in the ambient space;
  --   h        homotopy in n+1 variables;
  --   xt       start solution with its last component equal to zero,
  --            satisfies the homotopy h (upto tolerance);
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   xt       solution at the end of the path, tracked to the
  --            last component of xt to be equal to one;
  --   sol      standard representation of the solution.

-- TRACKING MANY PATHS WITHOUT OUTPUT :

  procedure Run ( n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xtsols : in out Standard_Complex_Solutions.Solution_List;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 );
  procedure Run ( n : in integer32;
                  h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                  xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 );
  procedure Run ( n : in integer32;
                  h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                  xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                  sols : in out QuadDobl_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution in xt,
  --   as defined by the homotopy h, without output,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   n        number of variables in the ambient space;
  --   h        homotopy in n+1 variables;
  --   xtsols   start solutions with their last component equal to zero,
  --            satisfies the homotopy h (upto tolerance);
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   xtsols   solutions at the end of the path, tracked to the
  --            last component of vectors in xtsols to be equal to one;
  --   sols     standard representation of the solutions.

  procedure Run ( file : in file_type; n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xtsols : in out Standard_Complex_Solutions.Solution_List;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 );
  procedure Run ( file : in file_type; n : in integer32;
                  h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                  xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 );
  procedure Run ( file : in file_type; n : in integer32;
                  h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                  xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                  sols : in out QuadDobl_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks paths starting at the solution in xtsols,
  --   as defined by the homotopy h, 
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics,
  --            if omitted, then there is no intermediate output;
  --   n        number of variables in the ambient space;
  --   h        homotopy in n+1 variables;
  --   xtsols   start solutions with their last component equal to zero,
  --            satisfies the homotopy h (upto tolerance);
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   xtsols   solutions at the end of the path, tracked to the
  --            last component of vectors in xtsols to be equal to one;
  --   sols     standard representation of the solutions.

-- MULTITASKED TRACKING OF MANY PATHS WITHOUT OUTPUT :

  procedure Multitasked_Run
              ( nv,nt : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out Standard_Complex_Solutions.Solution_List;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Multitasked_Run
              ( nv,nt : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Multitasked_Run
              ( nv,nt : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution in xt,
  --   as defined by the homotopy h, without output,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   nv       number of variables in the ambient space;
  --   nt       number of tasks for the multitasked path trackers;
  --   h        homotopy in nv+1 variables;
  --   xtsols   start solutions with their last component equal to zero,
  --            satisfies the homotopy h (upto tolerance).

  -- ON RETURN :
  --   xtsols   solutions at the end of the path, tracked to the
  --            last component of vectors in xtsols to be equal to one;
  --   sols     standard representation of the solutions.

-- MULTITASKED TRACKING OF MANY PATHS WITH OUTPUT TO FILE :

  procedure Multitasked_Run
              ( file : in file_type; nv,nt : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out Standard_Complex_Solutions.Solution_List;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Multitasked_Run
              ( file : in file_type; nv,nt : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Multitasked_Run
              ( file : in file_type; nv,nt : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Tracks paths starting at the solution in xtsols,
  --   as defined by the homotopy h, 
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics,
  --            if omitted, then there is no intermediate output;
  --   nv       number of variables in the ambient space;
  --   nt       number of tasks for the multitasked path trackers;
  --   h        homotopy in nv+1 variables;
  --   xtsols   start solutions with their last component equal to zero,
  --            satisfies the homotopy h (upto tolerance).

  -- ON RETURN :
  --   xtsols   solutions at the end of the path, tracked to the
  --            last component of vectors in xtsols to be equal to one;
  --   sols     standard representation of the solutions.

end Wrapped_Path_Trackers;
