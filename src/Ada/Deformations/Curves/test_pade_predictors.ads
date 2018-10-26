with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Vectors;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
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

  procedure Standard_Step_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Dense_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector );
  procedure DoblDobl_Step_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Dense_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector );
  procedure QuadDobl_Step_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Dense_Series_Vectors.Vector;
                pv : QuadDobl_Pade_Approximants.Pade_Vector );

  -- DESCRIPTION :
  --   Tests the determination of the step size interactively
  --   to gain computational experience with its effectiveness,
  --   in double, double double, or quad double precision.
  --   Compares the series predictor with the Pade predictor.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   pv       a vector of Pade approximants.


  -- DESCRIPTION :
  --   Tests the determination of the step size interactively
  --   to gain computational experience with its effectiveness,
  --   in double double precision.
  --   Compares the series predictor with the Pade predictor.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   pv       a vector of Pade approximants.

  procedure Standard_Test_Pade_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Solutions.Solution;
                nit : in integer32;
                pv : out Standard_Pade_Approximants.Pade_Vector );
  procedure DoblDobl_Test_Pade_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Solutions.Solution;
                nit : in integer32;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector );
  procedure QuadDobl_Test_Pade_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Solutions.Solution;
                nit : in integer32;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector );

  -- DESCRIPTION :
  --   Given a homotopy hom and a start solution sol,
  --   prompts the user for the degrees of the numerator and denominator
  --   of the Pade approximants, followed by tests on the Pade predictor,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   hom      a homotopy suitable for power series methods;
  --   sol      a solution for t = 0;
  --   nit      number of iteration of the Newton power method.

  -- ON RETURN :
  --   pv       the vector of Pade approximant.

end Test_Pade_Predictors;
