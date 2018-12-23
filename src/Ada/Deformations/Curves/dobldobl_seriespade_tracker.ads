with DoblDobl_Complex_Poly_Systems;      use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;         use DoblDobl_Complex_Solutions;
with Homotopy_Continuation_Parameters;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Pade_Approximants;

package DoblDobl_SeriesPade_Tracker is

-- DESCRIPTION :
--   The package implements a path tracker with a next() method.
--   The path tracker applies a Pade predictor computed via Newton's
--   method on power series, in double double precision.

-- CONSTRUCTORS :

  procedure Init ( pars : in Homotopy_Continuation_Parameters.Parameters );

  -- DESCRIPTION :
  --   Stores the values for the homotopy continuation parameters in pars.

  procedure Init ( p,q : in Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Initializes the homotopy with target system p and start system q,
  --   using the gamma, initialized with the Init(pars) procedure.

  -- REQUIRED :
  --   The homotopy continuation parameters have been initialized
  --   with the previous Init procedure.

  procedure Init ( s : in Link_to_Solution );

  -- DESCRIPTION :
  --   Stores s as the current solution for tracking the next path.

-- PREDICTOR-CORRECTOR STAGE :

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

  function Get_Current_Series_Vector
    return DoblDobl_Complex_Series_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the current series vector.

  function Get_Current_Pade_Vector
    return DoblDobl_Pade_Approximants.Link_to_Pade_Vector;

  -- DESCRIPTION :
  --   Returns the current vector of Pade approximants.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the homotopy data and resets all data.

end DoblDobl_SeriesPade_Tracker;
