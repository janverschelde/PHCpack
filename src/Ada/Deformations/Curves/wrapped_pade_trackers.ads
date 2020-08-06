with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;

package Wrapped_Pade_Trackers is

-- DESCRIPTION :
--   Wrappers are provided to the path trackers with Pade approximants
--   as predictors.  The wrappers instantiate the homotopy, given as
--   a polynomial system where the continuation parameter is the last
--   variable in the polynomials of the system.
--   The design of the package follows the wrapped_path_trackers,
--   used in the Littlewood-Richardson homotopies.

-- SETUP OF PARAMETERS :

  procedure Set_Parameters;
  procedure Set_Parameters ( file : in file_type );

  -- DESCRIPTION :
  --   Sets the default values for the homotopy continuation parameters
  --   and then allows the user to interactively tune the parameters.
  --   The settings of the parameter values are saved and written to
  --   file if a file is provided.

-- TRACKING ONE PATH WITH OUTPUT TO FILE :

  procedure Run ( file : in file_type; n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution in xt,
  --   as defined by the homotopy h, with intermediate output.

  -- REQUIRED : Set_Parameters() was executed.

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

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the data for the settings of the parameter.

end Wrapped_Pade_Trackers;
