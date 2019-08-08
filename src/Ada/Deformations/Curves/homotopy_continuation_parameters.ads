with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Homotopy_Continuation_Parameters is

-- DESCRIPTION :
--   The homotopy continuation parameters collect the settings of the
--   parameters and tolerances for the path trackers with Pade approximants.

  type Parameters is record
    alpha : double_float;    -- tolerance on the residual of the prediction
    pbeta : double_float;    -- multiplication factor of the pole radius
    cbeta : double_float;    -- multiplication factor for the curvature
    gamma : Complex_Number;  -- the gamma constant in the homotopy
    tolres : double_float;   -- tolerance on the residual in the corrector
    epsilon : double_float;  -- tolerance on zero series coefficients
    numdeg : natural32;      -- degree of numerator of Pade approximant
    dendeg : natural32;      -- degree of denominator of Pade approximant
    maxsize : double_float;  -- the maximum step size
    minsize : double_float;  -- the minimum step size
    corsteps : natural32;    -- maximum number of corrector steps
    maxsteps : natural32;    -- maximum number of steps on a path
  end record;

  type Link_to_Parameters is access Parameters;

  function Default_Values return Parameters;

  -- DESCRIPTION :
  --   Returns default values for the parameters.

  procedure Clear ( pars : in out Link_to_Parameters );

  -- DESCRIPTION :
  --   Deallocates the space for the parameters.

-- STORAGE of an instance of the parameters :

  procedure Construct ( pars : in Parameters );

  -- DESCRIPTION :
  --   Makes a pointer to pars which points to the parameters in pars.

  function Retrieve return Link_to_Parameters;

  -- DESCRIPTION :
  --   Returns the pointer stored with the Construct procedure.

  procedure Destruct;

  -- DESCRIPTION :
  --   Deallocates the space allocated with the Construct.

end Homotopy_Continuation_Parameters;
