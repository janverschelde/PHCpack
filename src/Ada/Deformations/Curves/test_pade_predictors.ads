with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Vectors;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_Systems;
with Standard_Pade_Approximants;
with DoblDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants;

package Test_Pade_Predictors is

-- DESCRIPTION :
--   This package provides procedures to test a Pade predictor.
--   A Pade predictor computes the Pade approximant of a solution
--   to a polynomial homotopy and applies the Pade approximant
--   to predict the next solution on the path.

  procedure Standard_Check_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Dense_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                step : in double_float );
  procedure DoblDobl_Check_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Dense_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                step : in double_double );
  procedure QuadDobl_Check_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Dense_Series_Vectors.Vector;
                pv : in QuadDobl_Pade_Approximants.Pade_Vector;
                step : in quad_double );

  -- DESCRIPTION :
  --   Checks the quality of the predicted step in double,
  --   double double, or quad double precision precision.
  --   Prints the evaluated Pade approximant and the norm of the
  --   evaluated predicted solution to screen.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   pv       vector of Pade approximants;
  --   step     a step size.

end Test_Pade_Predictors;
