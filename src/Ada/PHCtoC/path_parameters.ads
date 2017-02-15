with text_io;                         use text_io;
with Standard_Integer_Numbers;        use Standard_Integer_Numbers;
with Standard_Floating_Numbers;       use Standard_Floating_Numbers;

package Path_Parameters is

-- DESCRIPTION :
--   This package provides default values, tuning and output procedures
--   for the numerical parameters and tolerances of the path trackers
--   with algorithmic differentiation in the Path library.

-- DEFAULT VALUES :

  N_PREDICTOR : integer32 := 4;
  STEP_INCREASE : double_float := 1.25;
  STEP_DECREASE : double_float := 0.7;
  MAX_DELTA_T : double_float := 1.0E-1;
  MAX_DELTA_T_END : double_float := 1.0E-2;
  MIN_DELTA_T : double_float := 1.0E-7;
  ERR_MAX_RES : double_float := 1.0E-6;
  ERR_MAX_DELTA_X : double_float := 1.0E-6;
  ERR_MAX_FIRST_DELTA_X : double_float := 1.0E-2;

  D_MAX_STEP  : integer32 := 1000;
  DD_MAX_STEP : integer32 := 2000;
  QD_MAX_STEP : integer32 := 3000;

  D_MAX_IT : integer32 := 3;
  DD_MAX_IT : integer32 := 4;
  QD_MAX_IT : integer32 := 5;
  D_ERR_MIN_ROUND_OFF : double_float := 1.0E-9;
  DD_ERR_MIN_ROUND_OFF : double_float := 1.0E-14;
  QD_ERR_MIN_ROUND_OFF : double_float := 1.0E-26;

  D_MAX_IT_REFINE : integer32 := 3;
  DD_MAX_IT_REFINE : integer32 := 4;
  QD_MAX_IT_REFINE : integer32 := 5;
  D_ERR_MIN_ROUND_OFF_REFINE : double_float := 1.0E-11;
  DD_ERR_MIN_ROUND_OFF_REFINE : double_float := 1.0E-22;
  QD_ERR_MIN_ROUND_OFF_REFINE : double_float := 1.0E-40;

  type Parameters is record -- holds the 14 Path parameters
    max_step : integer32;           --  1. Maximum number of steps
    n_predictor : integer32;        --  2. Number of points in the predictor
    step_increase : double_float;   --  3. Increase factor on the step size
    step_decrease : double_float;   --  4. Decrease factor on the step size
    max_delta_t : double_float;     --  5. Maximal step size along a path
    max_delta_t_end : double_float; --  6. Maximal step size at end of a path
    min_delta_t : double_float;     --  7. Minimum step size along a path 
    err_max_res : double_float;     --  8. Tolerance on the residual
    err_max_delta_x : double_float; --  9. Tolerance on the corrector update
    err_max_first_delta_x : double_float;
                                    -- 10. Tolerance on the 1st corrector update
    max_it : integer32;             -- 11. Maximum number of Newton iterations
    err_min_round_off : double_float;
                                    -- 12. Tolerance for Newton's corrector
    max_it_refine : integer32;      -- 13. Maximum number of Newton refine steps
    err_min_round_off_refine : double_float;
                                    -- 14. Tolerance for Newton's refinement
  end record;

  function Default_Parameters ( precision : integer32 ) return Parameters;

  -- DESCRIPTION :
  --   Returns the default values for the parameters,
  --   for precision equal to 16, 32, or 64.

  procedure Write ( pars : in Parameters );
  procedure Write ( file : file_type; pars : in Parameters );

  -- DESCRIPTION :
  --   Writes the parameters to standard output or to file.

  procedure Set_Value ( pars : in out Parameters;
                        idx : in integer32; val : in double_float );

  -- DESCRIPTION :
  --   Given the index idx of the parameter, sets the corresponding
  --   value of the parameter to the value of val.

  procedure Tune ( pars : in out Parameters );

  -- DESCRIPTION :
  --   Interactive tuning of the parameters.

end Path_Parameters;
