with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_Systems;

package Series_and_Trackers is

-- DESCRIPTION :
--   Path trackers with Newton power series predictors are provided
--   in standard double, double double, or quad double precision.

  procedure Correct
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                t : in double_float; nit : in natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies a couple of Newton iterations to correct the solution sol
  --   of the system at t, using nit steps.

  procedure Track_one_Path
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                verbose : in boolean := false );
  procedure Track_one_Path
              ( h : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution;
                verbose : in boolean := false );
  procedure Track_one_Path
              ( h : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out Quaddobl_Complex_Solutions.Solution;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom,
  --   in standard double, double double, or quad double precision.

end Series_and_Trackers;
