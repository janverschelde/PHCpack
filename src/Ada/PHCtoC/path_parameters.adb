with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Numbers_io;

package body Path_Parameters is

  function Default_Parameters ( precision : integer32 ) return Parameters is

    res : Parameters;

  begin
    res.n_predictor := N_PREDICTOR;
    res.step_increase := STEP_INCREASE;
    res.step_decrease := STEP_DECREASE;
    res.max_delta_t := MAX_DELTA_T;
    res.max_delta_t_end := MAX_DELTA_T_END;
    res.min_delta_t := MIN_DELTA_T;
    res.err_max_res := ERR_MAX_RES;
    res.err_max_delta_x := ERR_MAX_DELTA_X;
    res.err_max_first_delta_x := ERR_MAX_FIRST_DELTA_X;
    if precision <= 16 then
      res.max_step := D_MAX_STEP;
      res.max_it := D_MAX_IT;
      res.max_it_refine := D_MAX_IT_REFINE;
      res.err_min_round_off := D_ERR_MIN_ROUND_OFF;
      res.err_min_round_off_refine := D_ERR_MIN_ROUND_OFF_REFINE;
    elsif precision <= 32 then
      res.max_step := DD_MAX_STEP;
      res.max_it := DD_MAX_IT;
      res.max_it_refine := DD_MAX_IT_REFINE;
      res.err_min_round_off := DD_ERR_MIN_ROUND_OFF;
      res.err_min_round_off_refine := DD_ERR_MIN_ROUND_OFF_REFINE;
    else
      res.max_step := QD_MAX_STEP;
      res.max_it := QD_MAX_IT;
      res.max_it_refine := QD_MAX_IT_REFINE;
      res.err_min_round_off := QD_ERR_MIN_ROUND_OFF;
      res.err_min_round_off_refine := QD_ERR_MIN_ROUND_OFF_REFINE;
    end if;
    return res;
  end Default_Parameters;

  procedure Write ( pars : in Parameters ) is
  begin
    Write(standard_output,pars);
  end Write;

  procedure Write ( file : file_type; pars : in Parameters ) is
  begin
    put(file," 1. Maximum number of steps                   : ");
    put(file,pars.max_step,1); new_line(file);
    put(file," 2. Number of points in the predictor         : ");
    put(file,pars.n_predictor,1); new_line(file);
    put(file," 3. Increase factor on the step size          :");
    put(file,pars.step_increase,2); new_line(file);
    put(file," 4. Decrease factor on the step size          :");
    put(file,pars.step_decrease,2); new_line(file);
    put(file," 5. Maximal step size along a path            :");
    put(file,pars.max_delta_t,2); new_line(file);
    put(file," 6. Maximal step size at the end of a path    :");
    put(file,pars.max_delta_t_end,2); new_line(file);
    put(file," 7. Minimum step size along a path            :");
    put(file,pars.min_delta_t,2); new_line(file);
    put(file," 8. Tolerance on the residual                 :");
    put(file,pars.err_max_res,2); new_line(file);
    put(file," 9. Tolerance on the corrector update         :");
    put(file,pars.err_max_delta_x,2); new_line(file);
    put(file,"10. Tolerance on the first corrector update   :");
    put(file,pars.err_max_first_delta_x,2); new_line(file);
    put(file,"11. Maximum number of Newton iterations       : ");
    put(file,pars.max_it,1); new_line(file);
    put(file,"12. Tolerance for Newton's corrector method   :");
    put(file,pars.err_min_round_off,2); new_line(file);
    put(file,"13. Maximum number of Newton refinement steps : ");
    put(file,pars.max_it_refine,1); new_line(file);
    put(file,"14. Tolerance for Newton's refinement method  :");
    put(file,pars.err_min_round_off_refine,2); new_line(file);
  end Write;

  procedure Set_Value ( pars : in out Parameters;
                        idx : in integer32; val : in double_float ) is
  begin
    case idx is
      when  1 => pars.max_step := integer32(val);
      when  2 => pars.n_predictor := integer32(val);
      when  3 => pars.step_increase := val;
      when  4 => pars.step_decrease := val;
      when  5 => pars.max_delta_t := val;
      when  6 => pars.max_delta_t_end := val;
      when  7 => pars.min_delta_t := val;
      when  8 => pars.err_max_res := val;
      when  9 => pars.err_max_delta_x := val;
      when 10 => pars.err_max_first_delta_x := val;
      when 11 => pars.max_it := integer32(val);
      when 12 => pars.err_min_round_off := val;
      when 13 => pars.max_it_refine := integer32(val);
      when 14 => pars.err_min_round_off_refine := val;
      when others => null;
    end case;
  end Set_Value;

  procedure Tune ( pars : in out Parameters ) is

    idx,natval : natural32 := 0;
    dblval : double_float := 0.0;

  begin
    loop
      new_line;
      put_line("The current values for the path parameters :");
      Write(pars);
      put("Type integer in 1..14 to change, 0 to exit : ");
      Numbers_io.Read_Natural(idx);
      exit when idx = 0;
      case idx is
        when 1 =>
          put("Give a new value for the maximum number of steps : ");
          Numbers_io.Read_Positive(integer(natval));
          pars.max_step := integer32(natval);
        when 2 =>
          put("Give a new value for the number of points in the predictor : ");
          Numbers_io.Read_Natural(natval);
          pars.n_predictor := integer32(natval);
        when 3 =>
          put("Give a new value for the increase factor on the step size : ");
          Numbers_io.Read_Positive_Float(dblval);
          pars.step_increase := dblval;
        when 4 =>
          put("Give a new value for the decrease factor on the step size : ");
          Numbers_io.Read_Positive_Float(dblval);
          pars.step_decrease := dblval;
        when 5 =>
          put("Give a new value for the maximal step long a path : ");
          Numbers_io.Read_Positive_Float(dblval);
          pars.max_delta_t := dblval;
        when 6 =>
          put("Give a new value for the maximal step at the end of a path :");
          Numbers_io.Read_Positive_Float(dblval);
          pars.max_delta_t_end := dblval;
        when 7 =>
          put("Give a new value for the minimum step size along a path : ");
          Numbers_io.Read_Positive_Float(dblval);
          pars.min_delta_t := dblval;
        when 8 =>
          put("Give a new value for the tolerance on the residual : ");
          Numbers_io.Read_Positive_Float(dblval);
          pars.err_max_res := dblval;
        when 9 =>
          put("Give a new value for the tolerance on the corrector update : ");
          Numbers_io.Read_Positive_Float(dblval);
          pars.err_max_delta_x := dblval;
        when 10 =>
          put("Give a new value for the tolerance on the first update : ");
          Numbers_io.Read_Positive_Float(dblval);
          pars.err_max_first_delta_x := dblval;
        when 11 =>
          put("Give a new value for the maximum number of corrector steps : ");
          Numbers_io.Read_Positive(integer(natval));
          pars.max_it := integer32(natval);
        when 12 =>
          put("Give a new value for the tolerance for the corrector method :");
          Numbers_io.Read_Positive_Float(dblval);
          pars.err_min_round_off := dblval;
        when 13 =>
          put("Give a new value for the maximum number of refinement steps : ");
          Numbers_io.Read_Positive(integer(natval));
          pars.max_it_refine := integer32(natval);
        when 14 =>
          put("Give a new value for the tolerance for the refinement :");
          Numbers_io.Read_Positive_Float(dblval);
          pars.err_min_round_off_refine := dblval;
        when others => null;
      end case;
    end loop;
  end Tune;

end Path_Parameters;
