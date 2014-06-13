with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Continuation_Parameters is

-- DESCRIPTION :
--   This package contains all important parameters and the default values.

-- DATA STRUCTURES FOR PARAMETERS :

  type Pred_Pars is record    -- contains the parameters for the predictor

    minstep,maxstep : double_float;  -- minimum and maximum step size
    expfac,redfac : double_float;    -- expansion and reduction factor
                                     --  for step length control
    success_steps : natural32;   -- number of successful steps before expansion
    predictor_type : natural32;  -- type of predictor used
    dist_target : double_float;  -- distance to target
    power : natural32;           -- power of t in (polyhedral) homotopy

  end record;

  type Corr_Pars is record    -- contains the parameters for the corrector

    epsrx,epsax,epsrf,epsaf : double_float;  
                              -- desired precisions for x and its residual f(x)
                              -- once relative (r) and once absolute (a)

    maxit,maxtot : natural32; -- maximum number of corrector iterations
                              -- for one step and for the whole path
  end record;

-- GLOBAL MONITOR :

  working_precision : natural32 := 16;

  condition : natural32 := 0;                                           --  1
            -- low = smooth paths, high = difficult path

  block_size : natural32 := 1;                                          --  2
            -- number of paths tracked simultaneously:
            -- sequential, parallel or a combination of both.

  max_steps : natural32 := 500*(condition+1);                           --  3
            -- maximum number of steps along a path

  start_end_game : double_float := 0.1;                                 --  4
            -- distance from target to start the end game

  endext_order : natural32 := 0;                                        --  5
            -- order of extrapolator in polyhedral end game

  max_reruns : natural32 := condition + 1;                              --  6
            -- maximum number of re-runs allowed

-- STEP CONTROL (PREDICTOR) :

  predictor_path_type : natural32 := 8;                                 --  7
  predictor_endg_type : natural32 := 8;                                 --  8
            -- predictor types:
            --   0 : secant for x, real for t;
            --   1 : secant for x, complex for t;
            --   2 : secant for x, geometric for t;
            --   3 : tangens for x, real for t;
            --   4 : tangens for x, complex for t;
            --   5 : tangens for x, geometric for t;
            --   6 : Hermite for x, real for t;
            --   7 : quadratic for x, real for t
            --   8 : cubic for x, real for t
            --   9 : quartic for x, real for t

  min_path_step_size : double_float
                     := 10.0**(-6 - integer(condition)/4);              --  9
  min_endg_step_size : double_float
                     := 10.0**(-8 - integer(condition)/4);              -- 10
            -- minimum step size along a path and at end of path

  max_path_step_size : double_float := 0.1/(double_float(condition+1));  -- 11
  max_endg_step_size : double_float := 0.01/(double_float(condition+1)); -- 12
            -- maximum step size along a path and at end of path

  reduction_path_factor : double_float := 0.70;                         -- 13
  reduction_endg_factor : double_float := 0.50;                         -- 14
            -- reduction factor for step size along a path and at end

  expansion_path_factor : double_float := 1.25;                         -- 15
  expansion_endg_factor : double_float := 1.10;                         -- 16
            -- expansion factor for step size

  success_path_steps : natural32 := 1 + 2*condition;                    -- 17
  success_endg_steps : natural32 := 1 + 4*condition;                    -- 18
            -- threshold on number of successful steps
            -- before expanding step size

  power_of_t : natural32 := 1;   -- power of t in the (polyhedral) homotopy

-- PATH CLOSENESS (CORRECTOR) :

  max_path_iter : natural32 := 3;                                       -- 19
  max_endg_iter : natural32 := 3;                                       -- 20
            -- maximum number of iterations for one corrector step
            -- along the path and at the end of the path

  relative_path_residual : double_float
                         := 10.0**(-9 - integer(condition));          -- 21
  relative_endg_residual : double_float
                         := 10.0**(-11 - integer(condition));         -- 22
  absolute_path_residual : double_float
                         := 10.0**(-9 - integer(condition));          -- 23
  absolute_endg_residual : double_float
                         := 10.0**(-11 - integer(condition));         -- 24
            -- desired precision for residuals along the path and at end
            --  |F(x)| < *_residual, once relative and once absolute

  relative_path_correction : double_float
                           := 10.0**(-9 - integer(condition));        -- 25
  relative_endg_correction : double_float
                           := 10.0**(-11 - integer(condition));       -- 26
  absolute_path_correction : double_float
                           := 10.0**(-9 - integer(condition));        -- 27
  absolute_endg_correction : double_float
                           := 10.0**(-11 - integer(condition));       -- 28
            -- desired precision for corrections along the path and at end
            --  |delta(x)| < *_correction, once relative and once absolute

-- SOLUTIONS (TOLERANCES) :

  tol_path_inverse_condition : double_float := 10.0**(-4);              -- 29
  tol_endg_inverse_condition : double_float := 10.0**(-12);             -- 30
            -- tolerance for inverse condition of jacobian to
            -- decide whether solution is singular or not

  tol_path_distance : double_float := 10.0**(-4);                       -- 31
  tol_endg_distance : double_float := 10.0**(-12);                      -- 32
            -- tolerance for two solutions x1, x2 to be clustered,
            -- when |x1(k) - x2(k)| < tol_*_distance, for all k

  tol_path_at_infinity : double_float := 10.0**8;                       -- 33
  tol_endg_at_infinity : double_float := 10.0**12;                      -- 34
            -- tolerance for a solution x to lie a infinity,
            -- when |x(k)| > tol_at_infinity, for a certain k

-- TUNING OF THE PARAMETERS :

  procedure Tune ( estimate_for_condition : in natural32 );
  procedure Tune ( estimate_for_condition,decimal_places : in natural32 );

  -- DESCRIPTION :
  --   Given an estimate for the condition of the homotopy, parameters
  --   will be set, according to the formulas in the defaults.
  --   The number of decimal places allows for multi-precision tracing.

  procedure Tune_Endgm_Pred ( extrapolation_order : in natural32 );

  -- DESCRIPTION :
  --   Determines the settings for the predictor in the end game,
  --   depending on the extrapolation order.

-- CREATING PARAMETER SETS :

  function Create_for_Path return Pred_Pars;
  function Create_End_Game return Pred_Pars;
  function Create_for_Path return Corr_Pars;
  function Create_End_Game return Corr_Pars;

  -- DESCRIPTION :
  --   Given the values for the parameters, sets of predictor and
  --   corrector parameters will be created.

end Continuation_Parameters;
