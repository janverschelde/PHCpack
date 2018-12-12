with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Homotopy_Continuation_Parameters;

package Standard_SeriesPade_Tracker is

-- DESCRIPTION :
--   The package implements a path tracker with a next() method.
--   The path tracker applies a Pade predictor computed via Newton's
--   method on power series.

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

-- SELECTORS :

  function Get_Parameters
    return Homotopy_Continuation_Parameters.Link_to_Parameters;

  -- DESCRIPTION :
  --   Returns the link to the current values of the parameters.
  --   This link can be used to change the values of the parameters.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the homotopy data and resets all data.

end Standard_SeriesPade_Tracker;
