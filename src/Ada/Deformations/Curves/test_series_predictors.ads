with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with PentDobl_Complex_Solutions;
with OctoDobl_Complex_Solutions;
with DecaDobl_Complex_Solutions;
with HexaDobl_Complex_Solutions;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;

package Test_Series_Predictors is

-- DESCRIPTION :
--   This package provides procedures to test a series predictor.
--   A series predictor computes the power series approximation of
--   a solution to a polynomial homotopy and applies the series
--   to predict the next solution on the path.

  procedure Standard_Check_Prediction
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Complex_Series_Vectors.Vector;
                step : in double_float );
  procedure DoblDobl_Check_Prediction
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Complex_Series_Vectors.Vector;
                step : in double_double );
  procedure QuadDobl_Check_Prediction
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Complex_Series_Vectors.Vector;
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
              ( nbeq,idxpar : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure TripDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out TripDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List );
  procedure PentDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out PentDobl_Complex_Solutions.Solution_List );
  procedure OctoDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out OctoDobl_Complex_Solutions.Solution_List );
  procedure DecaDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out DecaDobl_Complex_Solutions.Solution_List );
  procedure HexaDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out HexaDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Wraps the reading of a homotopy in double, double double,
  --   triple double, quad double, penta double, octo double,
  --   deca double, or hexa double precision.
  --
  -- ON RETURN :
  --   nbeq     the number of equations in the homotopy;
  --   idxpar   index of the natural parameter, which is zero
  --            in case the homotopy has an artificial parameter;
  --   sols     the start solutions.

end Test_Series_Predictors;
