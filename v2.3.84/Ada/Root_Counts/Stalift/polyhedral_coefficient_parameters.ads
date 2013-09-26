with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Polyhedral_Coefficient_Parameters is

-- DESCRIPTION :
--   This package contains the important numerical parameters 
--   for the polyhedral coefficient path trackers,
--   with their default values and tuning facilities.

-- PARAMETERS FOR PREDICTOR :

  min_infinity : double_float := -4.5;
    -- starting value for the s parameter in the exponential,
    -- while exp(-infinity) = 0, it suffices to take a number
    -- sufficiently negative to have a good starting value

  max_pred_step : double_float := 2.0;
    -- maximal predictor step for s, the independent continuation
    -- parameter, this value is linked to the value of min_infinity:
    -- the more negative min_infinity, the larger max_pred_step may be

  expansion_threshold : natural32 := 5;
    -- required number of consecutive successful corrector stages
    -- before the step size is expanded

  expansion_factor : double_float := 1.5;
    -- expansion factor to multiply the step size with
    -- in case of a successful correction

  reduction_factor : double_float := 0.5;
    -- reduction factor to multiply the step size with
   --  in case of a failed correction

-- PARAMETERS FOR CORRECTOR :

  max_corr_iter : natural32 := 4;
    -- maximal number of corrector iterations

  tol_root : double_float := 1.0E-9;
    -- tolerance on magnitude of correction term and residual
    -- while tracking a path

  max_nb_steps : natural32 := 500;
    -- maximal number of steps along a path

-- TUNING FACILITIES :

  procedure Write ( file : in file_type );

  -- DESCRIPTION :
  --   Writes the current values for the parameters to file.

  procedure Tune;

  -- DESCRIPTION :
  --   Interactive facility to tune the parameters.

end Polyhedral_Coefficient_Parameters;
