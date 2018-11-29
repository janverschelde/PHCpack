with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Homotopy_Continuation_Parameters is

-- DESCRIPTION :
--   The homotopy continuation parameters collect the settings of the
--   parameters and tolerances for the path trackers with Pade approximants.

  type Parameters is record
    alpha : double_float;    -- tolerance on the residual of the prediction
    sbeta : double_float;    -- multiplication factor of the series step
    pbeta : double_float;    -- multiplication factor of the pole radius
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

  function Default_Values return Parameters;

  -- DESCRIPTION :
  --   Returns default values for the parameters.

end Homotopy_Continuation_Parameters;
