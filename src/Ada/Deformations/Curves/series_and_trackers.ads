with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;  
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Homotopy_Continuation_Parameters;

package Series_and_Trackers is

-- DESCRIPTION :
--   Path trackers with Newton power series predictors are provided
--   in standard double, double double, or quad double precision.
--   The versions may be silent or verbose.

  procedure Correct
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                t : in double_float; nit : in natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float );
  procedure Correct
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in double_double; nit : in natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double );
  procedure Correct
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in quad_double; nit : in natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double );

  -- DESCRIPTION :
  --   Applies Newton's method to correct the solution, silent version,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   hom      the homotopy with series coefficients;
  --   t        value for th series parameter to evaluate the series in hom;
  --   nit      number of steps to do with Newton's method;
  --   sol      predicted value for the solution.

  -- ON RETURN :
  --   err      magnitude of the correction term;
  --   rco      estimate for the inverse condition number;
  --   res      magnitude of the residual.

  procedure Correct
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                t : in double_float; nit : in natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float;
                verbose : in boolean := false );
  procedure Correct
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in double_double; nit : in natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in boolean := false );
  procedure Correct
              ( file : in file_type;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in quad_double; nit : in natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies Newton's method to correct the solution, verbose version,
  --   in standard double or double double precision.

  -- ON ENTRY :
  --   file     for writing extra diagnostic output, if verbose;
  --   hom      the homotopy with series coefficients;
  --   t        value for th series parameter to evaluate the series in hom;
  --   nit      number of steps to do with Newton's method;
  --   sol      predicted value for the solution;
  --   verbose  to indicate that extra output is wanted.

  -- ON RETURN :
  --   err      magnitude of the correction term;
  --   rco      estimate for the inverse condition number;
  --   res      magnitude of the residual.

  procedure Set_Step
              ( t,step : in out double_float;
                maxstep,target : in double_float );

  -- DESCRIPTION :
  --   Sets the step size and the new value for t,
  --   taking into account the maximal step size maxstep
  --   and the value for the target for t.

  function Residual_Prediction
              ( hom : Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : Standard_Complex_Vectors.Vector;
                step : double_float ) return double_float;
  function Residual_Prediction
              ( hom : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : DoblDobl_Complex_Vectors.Vector;
                step : double_float ) return double_float;
  function Residual_Prediction
              ( hom : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : QuadDobl_Complex_Vectors.Vector;
                step : double_float ) return double_float;

  -- DESCRIPTION :
  --   Given a homotopy, a predicted solution, with a step,
  --   returns the norm of the evaluated homotopy at the solution.

  procedure Track_one_Path
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure Track_one_Path
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure Track_one_Path
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom,
  --   in standard double, double double, or quad double precision.
  --   This version remains silent and does not write any output.

  -- ON ENTRY :
  --   hom      a homotopy with series coefficients;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances.

  -- ON RETURN :
  --   sol      solution at the end of the paths.

  procedure Track_one_Path
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );
  procedure Track_one_Path
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );
  procedure Track_one_Path
              ( file : in file_type;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom,
  --   in standard double, double double, or quad double precision.
  --   This version is verbose and writes extra diagnostic output to file.

  -- ON ENTRY :
  --   file     for writing output during the computations;
  --   hom      a homotopy with series coefficients;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances.

  -- ON RETURN :
  --   sol      solution at the end of the paths.

  procedure Track_Many_Paths
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure Track_Many_Paths
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure Track_Many_Paths
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );

  -- DESCRIPTION :
  --   Tracks the paths starting at the solutions in the list sols,
  --   as defined by the homotopy in hom,
  --   in double, double double, or quad double precision.
  --   The procedures are silent and do not write any output.

  -- ON ENTRY :
  --   hom      a homotopy with series coefficients;
  --   sols     start solutions in the homotopy;
  --   pars     values of the parameters and tolerances.

  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );
  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );
  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks the paths starting at the solutions in the list sols,
  --   as defined by the homotopy in hom,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     for output to file;
  --   hom      a homotopy with series coefficients;
  --   sols     start solutions in the homotopy;
  --   pars     values of the parameters and tolerances;
  --   verbose  if true, then extra output is written to file.

  -- ON RETURN :
  --   sols     solutions at the end of the paths.

end Series_and_Trackers;
