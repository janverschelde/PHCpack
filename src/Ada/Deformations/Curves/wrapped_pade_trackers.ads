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

-- TRACKING ONE PATH WITH OR WITHOUT OUTPUT TO FILE :

  procedure Run ( idxpar : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 );
  procedure Run ( file : in file_type; idxpar : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  verbose : in boolean := false;
                  vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution in xt,
  --   as defined by the homotopy h, with intermediate output.

  -- REQUIRED : Set_Parameters() was executed.

  -- ON ENTRY :
  --   file     optional output file for diagnostics and statistics,
  --            if omitted, then there is no intermediate output;
  --   idxpar   the continuation parameter is one of the variables in h,
  --            idxpar is the index of the variable in h that plays
  --            the role of the continuation parameter;
  --   h        a natural parameter homotopy;
  --   xt       start solution which statisfies h if the value
  --            for the variable with index equal to idxpar is zero;
  --   verbose  if extra output is needed (only possible with file);
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   xt       solution at the end of the path, tracked to the
  --            last component of xt to be equal to one;
  --   sol      standard representation of the solution.

-- TRACKING MANY PATHS WITH OR WITHOUT OUTPUT TO FILE :

  procedure Run ( idxpar : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 );
  procedure Run ( file : in file_type; idxpar : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  verbose : in boolean := false;
                  vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks many paths starting at the solution in sols,
  --   as defined by the homotopy h, with intermediate output.

  -- REQUIRED : Set_Parameters() was executed.

  -- ON ENTRY :
  --   file     optional output file for diagnostics and statistics,
  --            if omitted, then there is no intermediate output;
  --   idxpar   the continuation parameter is one of the variables in h,
  --            idxpar is the index of the variable in h that plays
  --            the role of the continuation parameter;
  --   h        a natural parameter homotopy;
  --   xtsols   start solutions which satisfy h if the value for
  --            the variable with index equal to idxpar is zero;
  --   verbose  if extra output is needed (only possible with file);
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   xtsols   solutions at the end of the path, tracked to the
  --            variable with index equal to idxpar in h, set to one;
  --   sols     standard representation of the solutions.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the data for the settings of the parameter.

end Wrapped_Pade_Trackers;
