with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Vectors;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_Systems;

package Test_Series_Predictors is

-- DESCRIPTION :
--   This package provides procedures to test a series predictor.
--   A series predictor computes the power series approximation of
--   a solution to a polynomial homotopy and applies the series
--   to predict the next solution on the path.

  procedure Standard_Check_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Dense_Series_Vectors.Vector;
                step : in double_float );
  procedure DoblDobl_Check_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Dense_Series_Vectors.Vector;
                step : in double_double );
  procedure QuadDobl_Check_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Dense_Series_Vectors.Vector;
                step : in quad_double );

  -- DESCRIPTION :
  --   Checks the quality of the predicted step in standard double,
  --   double double, or quad double precision.
  --   Writes the predicted solution and the norm of its evaluation.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   step     a step size.

  procedure Standard_Homotopy_Reader
              ( nbeq : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Homotopy_Reader
              ( nbeq : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Homotopy_Reader
              ( nbeq : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Wraps the reading of a homotopy in double, double double,
  --   or quad double precision.
  --
  -- ON RETURN :
  --   nbeq     the number of equations in the homotopy;
  --   sols     the start solutions.

end Test_Series_Predictors;
