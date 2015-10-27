with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;

package QuadDobl_Quad_Parameters is

-- DESCRIPTION :
--   This package collects the numerical parameters to tune the path
--   trackers to look for singular solutions in a sweeping homotopy.
--   For coding ease, all parameters are expressed as quad doubles.
--   Since double double arithmetic tolerates worser condition numbers,
--   the tolerances may thus then also be adjusted for harder problems.

-- STEP CONTROL for PREDICTOR :

  max_step_size : quad_double := create(0.1);
  reduction_multiplier : quad_double := create(0.5);
  expansion_multiplier : quad_double := create(1.5);
  expansion_threshold : natural32 := 3;

-- SOLUTION TOLERANCES for CORRECTOR :

  increment_tolerance : quad_double := create(1.0E-8);
  residual_tolerance : quad_double := create(1.0E-8);
  max_corrector_steps : natural32 := 4;

-- TO DECIDE WHEN CRITICAL :

  determinant_tolerance : quad_double := create(1.0E-12);

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

end QuadDobl_Quad_Parameters;
