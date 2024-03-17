with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Poly_Systems;      use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;
with Homotopy_Continuation_Parameters;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Pade_Approximants;

package QuadDobl_SeriesPade_Tracker is

-- DESCRIPTION :
--   The package implements a path tracker with a next() method.
--   The path tracker applies a Pade predictor computed via Newton's
--   method on power series, in quad double precision.

-- CONSTRUCTORS :

  procedure Init ( pars : in Homotopy_Continuation_Parameters.Parameters );

  -- DESCRIPTION :
  --   Stores the values for the homotopy continuation parameters in pars.

  procedure Init ( p,q : in Link_to_Poly_Sys; homogeneous : in boolean );

  -- DESCRIPTION :
  --   Initializes the homotopy with target system p and start system q,
  --   using the gamma, initialized with the Init(pars) procedure.
  --   If homogeneous, then 1-homogenization is applied,
  --   otherwise the tracking happens in the original affine coordinates.

  -- REQUIRED :
  --   The homotopy continuation parameters have been initialized
  --   with the previous Init procedure.

  procedure Init ( h : in Link_to_Poly_Sys; idx : in integer32 );

  -- DESCRIPTION :
  --   Initializes the tracker with the natural parameter homotopy h,
  --   where the index of the continuation parameter equals idx,
  --   the index in the list of variables of the homotopy.

  procedure Init ( s : in Link_to_Solution );

  -- DESCRIPTION :
  --   Stores s as the current solution for tracking the next path.

-- PREDICTOR-CORRECTOR STAGE :

  procedure Predictor_Feedback_Loop
              ( fail : out boolean; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs the predictor feedback loop with intermediate output if verbose.
  --   Each iteration of the loop computes the predictor residual and
  --   the step size in halved if the predictor residual is too large.
  --   If the predictor residual is small enough, then fail is false on return,
  --   otherwise, if the step size is smaller than the minimum step size,
  --   then fail will be false on return.

  procedure Predict ( fail : out boolean; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Starting at the current solution, computes the next power series,
  --   the next vector of Pade approximants, sets the step size, and
  --   sets the current solution to the predicted solution.
  --   If verbose, then extra output is written to screen,
  --   otherwise, the procedure remains silent.
  --   The fail flag is true on return if the step size dropped below
  --   the minimum step size set in the homotopy continuation parameters.

  procedure Correct ( fail : out boolean; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies Newton's method to correct the current solution.
  --   If verbose, then extra output is written to screen,
  --   otherwise, the procedure remains silent.
  --   On return is the fail flag, which is true if the required
  --   accuracy is not met within the allowed number of corrector steps,
  --   and/or when Newton's method diverges.

  procedure Predict_and_Correct
              ( fail : out boolean; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies the Predict and Correct stages after each other.
  --   If verbose, then extra output is written to screen,
  --   otherwise, the procedure remains silent.
  --   If fail is true on return, then the tracker failed to meet
  --   the requested accuracy requirements.

-- SELECTORS :

  function Get_Parameters
    return Homotopy_Continuation_Parameters.Link_to_Parameters;

  -- DESCRIPTION :
  --   Returns the link to the current values of the parameters.
  --   This link can be used to change the values of the parameters.

  function Get_Current_Solution return Link_to_Solution;

  -- DESCRIPTION :
  --   Returns the current solution.

  function Get_Predicted_Solution return Link_to_Solution;

  -- DESCRIPTION :
  --   Returns the predicted solution.

  function Get_Current_Series_Vector
    return QuadDobl_Complex_Series_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the current series vector.

  function Get_Current_Pade_Vector
    return QuadDobl_Pade_Approximants.Link_to_Pade_Vector;

  -- DESCRIPTION :
  --   Returns the current vector of Pade approximants.

  function Get_Current_Poles return QuadDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the roots of all the denominators of the Pade vector.
  --   The i-th component in the vecvec on return contains the poles
  --   for the i-th Pade approximant.

  function Get_Current_Pole_Radius return quad_double;

  -- DESCRIPTION :
  --   Returns the smallest pole radius, computed by the predictor.

  function Get_Current_Closest_Pole
             return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the closest pole, computed by the predictor.

  function Get_Current_Series_Step return double_float;

  -- DESCRIPTION :
  --   Returns the current value of the series step.

  function Get_Current_Pole_Step return double_float;

  -- DESCRIPTION :
  --   Returns the current value of the pole step.

  function Get_Current_Estimated_Distance return quad_double;

  -- DESCRIPTION :
  --   Returns the estimate distance to the closest solution.

  function Get_Current_Hessian_Step return double_float;

  -- DESCRIPTION :
  --   Returns the current value of the Hessian step.

  function Get_Current_Step_Size return double_float;

  -- DESCRIPTION :
  --   Returns the current value of the step size.

  function Get_Current_t_Value return double_float;

  -- DESCRIPTION :
  --   Returns the current value of the continuation parameter.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the homotopy data and resets all data.

end QuadDobl_SeriesPade_Tracker;
