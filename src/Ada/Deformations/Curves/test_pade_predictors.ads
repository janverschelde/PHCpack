with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
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
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                step : in double_float );
  procedure DoblDobl_Check_Prediction
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                step : in double_double );
  procedure QuadDobl_Check_Prediction
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
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
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Complex_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector );
  procedure DoblDobl_Step_Prediction
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Complex_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector );
  procedure QuadDobl_Step_Prediction
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Complex_Series_Vectors.Vector;
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

  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                solv : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector;
                pv : out Standard_Pade_Approximants.Pade_Vector;
                poles : out Standard_Complex_VecVecs.VecVec;
                fpr : out double_float; verbose : in boolean := false );
  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                solv : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector;
                poles : out DoblDobl_Complex_VecVecs.VecVec;
                fpr : out double_double; verbose : in boolean := false );
  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                solv : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector;
                poles : out QuadDobl_Complex_VecVecs.VecVec;
                fpr : out quad_double; verbose : in boolean := false );

  -- DESCRIPTION :
  --   The forward pole radius is the smallest radius of all poles
  --   of the vector of Pade approximants originating at the solution.

  -- REQUIRED :
  --   For double, double double, or quad double precision,
  --   the respective packages Standard_Homotopy, DoblDobl_Homotopy, 
  --   or QuadDobl_Homotopy have been initialized properly.

  -- ON ENTRY:
  --   neq      number of equations in the homotopy;
  --   degnum   degree of the numerator;
  --   degden   degree of the denominator;
  --   sol      the start solution for the homotopy;
  --   verbose  if true, then additional output is written to screen.

  -- ON RETURN :
  --   srv      series computed starting at the solution sol;
  --   eva      the evaluated solution series;
  --   pv       the vector of Pade approximants;
  --   poles    poles of the vector of Pade approximants;
  --   fpr      the smallest radius of all poles.

  function Forward_Pole_Radius
              ( neq,degnum,degden : integer32;
                solv : Standard_Complex_Vectors.Vector;
                verbose : boolean := false ) return double_float;
  function Forward_Pole_Radius
              ( neq,degnum,degden : integer32;
                solv : DoblDobl_Complex_Vectors.Vector;
                verbose : boolean := false ) return double_double;
  function Forward_Pole_Radius
              ( neq,degnum,degden : integer32;
                solv : QuadDobl_Complex_Vectors.Vector;
                verbose : boolean := false ) return quad_double;

  -- DESCRIPTION :
  --   This function wraps the corresponding procedure with the same name.
  --   It returns the smallest radius of all poles of the vector of Pade
  --   approximants computed at the given solution vector.

  -- REQUIRED :
  --   For double, double double, or quad double precision,
  --   the respective packages Standard_Homotopy, DoblDobl_Homotopy, 
  --   or QuadDobl_Homotopy have been initialized properly.

  -- ON ENTRY:
  --   neq      number of equations in the homotopy;
  --   degnum   degree of the numerator;
  --   degden   degree of the denominator;
  --   sol      the start solution for the homotopy;
  --   verbose  if true, then additional output is written to screen.

  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List;
                radii : out Standard_Floating_Vectors.Vector;
                fpr : out double_float; verbose : in boolean := false );
  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                radii : out Double_Double_Vectors.Vector;
                fpr : out double_double; verbose : in boolean := false );
  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                radii : out Quad_Double_Vectors.Vector;
                fpr : out quad_double; verbose : in boolean := false );

  -- DECRIPTION :
  --   For a given list of solutions, computes the forward pole radius
  --   for each solution and returns in fpr the smallest value.

  -- REQUIRED :
  --   For double, double double, or quad double precision,
  --   the respective packages Standard_Homotopy, DoblDobl_Homotopy, 
  --   or QuadDobl_Homotopy have been initialized properly.

  -- ON ENTRY:
  --   neq      number of equations in the homotopy;
  --   degnum   degree of the numerator;
  --   degden   degree of the denominator;
  --   sols     the start solutions for the homotopy;
  --   verbose  if true, then additional output is written to screen.

  -- ON RETURN :
  --   radii    an array of range 1..Length_Of(sols),
  --            radii(k) contains the forward pole radius for the
  --            k-th solution in sols;
  --   fpr      the smallest value over all radii.

  procedure Standard_Test_Pade_Prediction
              ( neq : in integer32;
                sol : in Standard_Complex_Solutions.Solution;
                pv : out Standard_Pade_Approximants.Pade_Vector );
  procedure DoblDobl_Test_Pade_Prediction
              ( neq : in integer32;
                sol : in DoblDobl_Complex_Solutions.Solution;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector );
  procedure QuadDobl_Test_Pade_Prediction
              ( neq : in integer32;
                sol : in QuadDobl_Complex_Solutions.Solution;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector );

  -- DESCRIPTION :
  --   Given a homotopy hom and a start solution sol,
  --   prompts the user for the degrees of the numerator and denominator
  --   of the Pade approximants, followed by tests on the Pade predictor,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   neq      number of equations in the homotopy;
  --   sol      a solution for t = 0.

  -- ON RETURN :
  --   pv       the vector of Pade approximant.

end Test_Pade_Predictors;
