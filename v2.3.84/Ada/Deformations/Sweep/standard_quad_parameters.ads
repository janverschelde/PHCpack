with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Standard_Quad_Parameters is

-- DESCRIPTION :
--   This package collects the numerical parameters to tune the path
--   trackers to look for singular solutions in a sweeping homotopy.

-- STEP CONTROL for PREDICTOR :

  max_step_size : double_float := 0.1;
  reduction_multiplier : double_float := 0.5;
  expansion_multiplier : double_float := 1.5;
  expansion_threshold : natural32 := 3;

-- SOLUTION TOLERANCES for CORRECTOR :

  increment_tolerance : double_float := 1.0E-8;
  residual_tolerance : double_float := 1.0E-8;
  max_corrector_steps : natural32 := 4;

-- TO DECIDE WHEN CRITICAL :

  determinant_tolerance : double_float := 1.0E-12;

-- GLOBAL STOP CRITERIUM :

  max_predictor_steps : natural32 := 500;

-- RESETTING to DEFAULTS :

  procedure Reset;

  -- DESCRIPTION :
  --   Resets the current values to the default values of above.

-- REPORTING and TUNING :

  procedure Show;
  procedure Show ( file : in file_type );

  -- DESCRIPTION :
  --   Writes the current settings of the parameters to screen or file.

  procedure Tune;

  -- DESCRIPTION :
  --   Interactive routine for user to change parameter settings.

end Standard_Quad_Parameters;
