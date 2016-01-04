package body Continuation_Parameters is

  procedure Tune ( estimate_for_condition : in natural32 ) is
  begin
    if working_precision > 16 then
      Tune(estimate_for_condition,working_precision);
    else
      condition  := estimate_for_condition;
      block_size := 1;
      max_steps  := 500*(condition + 1);
      max_reruns := condition + 1;
      predictor_path_type := 8; -- 0; cubic as default
      min_path_step_size := 10.0**(-6 - integer(condition)/4);
      max_path_step_size := 0.1/(double_float(condition+1));
      success_path_steps := 1 + condition;
      reduction_path_factor := 0.70 - 0.10*double_float(condition);
      reduction_endg_factor := 0.50 - 0.10*double_float(condition);
      if condition >= 4 then
        relative_path_residual   := 10.0**( -9 - integer(condition) mod 4);
        absolute_path_residual   := 10.0**( -9 - integer(condition) mod 4);
        relative_path_correction := 10.0**( -9 - integer(condition) mod 4);
        absolute_path_correction := 10.0**( -9 - integer(condition) mod 4);
        relative_endg_residual   := 10.0**(-11 - integer(condition) mod 4);
        absolute_endg_residual   := 10.0**(-11 - integer(condition) mod 4);
        relative_endg_correction := 10.0**(-11 - integer(condition) mod 4);
        absolute_endg_correction := 10.0**(-11 - integer(condition) mod 4);
      end if;
      Tune_Endgm_Pred(endext_order);
    end if;
  end Tune;

  procedure Tune ( estimate_for_condition,decimal_places : in natural32 ) is

    factor,dblfac,pathtol,endgtol : double_float;

  begin
    working_precision := decimal_places;
    tol_path_inverse_condition := 10.0**(-integer(decimal_places)/4);
    tol_endg_inverse_condition := 10.0**(-integer(decimal_places)/4-8);
    factor := tol_path_inverse_condition;
    dblfac := factor*factor;
    condition  := estimate_for_condition;
    block_size := 1;
    max_steps  := 500*(condition + 1);
    max_reruns := condition + 1;
    predictor_path_type := 8; -- 0; cubic as default
    min_path_step_size := 10.0**(-6 - integer(condition)/4);
    max_path_step_size := 0.1/(double_float(condition+1));
    if decimal_places >= 32 then -- first guess, needs better care !
      max_path_step_size := 0.001/(double_float(condition+1));
      max_steps := 50000*(condition+1);
    end if;
    success_path_steps := 1 + condition;
    reduction_path_factor := 0.70 - 0.10*double_float(condition);
    reduction_endg_factor := 0.50 - 0.10*double_float(condition);
    pathtol := 10.0**(-integer(decimal_places)-4*integer(condition))/dblfac;
    endgtol := 10.0**(-integer(decimal_places)-4*integer(condition))/factor;
    relative_path_residual   := pathtol;
    absolute_path_residual   := pathtol;
    relative_path_correction := pathtol;
    absolute_path_correction := pathtol;
    relative_endg_residual   := endgtol;
    absolute_endg_residual   := endgtol;
    relative_endg_correction := endgtol;
    absolute_endg_correction := endgtol;
    Tune_Endgm_Pred(endext_order);
  end Tune;

  procedure Tune_Endgm_Pred ( extrapolation_order : in natural32 ) is
  begin
    if extrapolation_order = 0 then
      min_endg_step_size := 10.0**(-8 - integer(condition)/4);
      max_endg_step_size := 0.01/(double_float(condition+1));
      success_endg_steps := 3 + condition;
    else
     -- if predictor_path_type <= 2
     --  then 
      predictor_endg_type := 2; -- no tangent predictor
     --  else predictor_endg_type := 5;
     -- end if;
      min_endg_step_size := min_path_step_size;
      max_endg_step_size := 0.5;
      reduction_endg_factor := 0.5;
      expansion_endg_factor := 1.7;
      success_endg_steps := 3 + condition;
    end if;
  end Tune_Endgm_Pred;

  function Create_for_Path return Pred_Pars is

    res : Pred_Pars;

  begin
    res.maxstep := max_path_step_size;
    res.minstep := min_path_step_size;
    res.expfac  := expansion_path_factor;
    res.redfac  := reduction_path_factor;
    res.success_steps := success_path_steps;
    res.predictor_type := predictor_path_type;
    res.dist_target := start_end_game;
    res.power := power_of_t;
    return res;
  end Create_for_Path;

  function Create_End_Game return Pred_Pars is

    res : Pred_Pars;

  begin
    res.maxstep := max_endg_step_size;
    res.minstep := min_endg_step_size;
    res.expfac  := expansion_endg_factor;
    res.redfac  := reduction_endg_factor;
    res.success_steps := success_endg_steps;
    res.predictor_type := predictor_endg_type;
    res.dist_target := 0.0;
    res.power := power_of_t;
    return res;
  end Create_End_Game;

  function Create_for_Path return Corr_Pars is

    res : Corr_Pars;

  begin
    res.epsrx  := relative_path_correction;
    res.epsax  := absolute_path_correction;
    res.epsrf  := relative_path_residual;
    res.epsaf  := absolute_path_residual;
    res.maxit  := max_path_iter;
    res.maxtot := max_steps*res.maxit;
    return res;
  end Create_for_Path;

  function Create_End_Game return Corr_Pars is

    res : Corr_Pars;

  begin
    res.epsrx  := relative_endg_correction;
    res.epsax  := absolute_endg_correction;
    res.epsrf  := relative_endg_residual;
    res.epsaf  := absolute_endg_residual;
    res.maxit  := max_endg_iter;
    res.maxtot := max_steps*res.maxit;
    return res;
  end Create_End_Game;

end Continuation_Parameters;
