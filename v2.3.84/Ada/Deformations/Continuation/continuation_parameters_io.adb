with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Numbers_io;                         use Numbers_io;
with Characters_and_Numbers;
with Continuation_Parameters;            use Continuation_Parameters;

package body Continuation_Parameters_io is

-- AUXILIARIES :

  function Predictor_Banner ( pred_type : natural32 ) return string is

  -- DESCRIPTION :
  --   Returns banner for the corresponding predictor type.

    s1s : constant string := "( x:Sec,t:Rea )";
    s1c : constant string := "( x:Sec,t:Com )";
    s1g : constant string := "( x:Sec,t:Geo )";
    t1s : constant string := "( x:Tan,t:Rea )";
    t1c : constant string := "( x:Tan,t:Com )";
    t1g : constant string := "( x:Tan,t:Geo )";
    h3s : constant string := "( x:Her,t:Rea )";
    q2s : constant string := "( x:Qu2,t:Rea )";
    c3s : constant string := "( x:Cub,t:Rea )";
    uuu : constant string := " no predictor  ";

  begin
    case pred_type is
      when 0 => return s1s;
      when 1 => return s1c;
      when 2 => return s1g;
      when 3 => return t1s;
      when 4 => return t1c;
      when 5 => return t1g;
      when 6 => return h3s;
      when 7 => return q2s;
      when 8 => return c3s;
      when others => return uuu;
    end case;
  end Predictor_Banner;

-- INFORMATION BALLOONS :

  procedure Condition_Info is

    i : array(1..4) of string(1..65);

  begin
    i(1):="  This global parameter  can  be  used  to  tune  the  parameters";
    i(2):="according to the general expectancy of how difficult the solution";
    i(3):="paths will be.  Low values of this parameter lead to a loose path";
    i(4):="tracking, whereas higher values tighten up the path-following.   ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Condition_Info;

  procedure Number_of_Paths_Info is

    i : array(1..9) of string(1..65);

  begin
    i(1):="  In simultaneous path-following, the same discretization for the";
    i(2):="values  of  the  continuation  parameter  is used for a number of";
    i(3):="paths.   Clustering  of  solution   paths   is   checked   during";
    i(4):="continuation   for   the   number  of  paths  that  are  followed";
    i(5):="simultaneously.                                                  ";
    i(6):="  In sequential path-following, i.e., when this parameter  equals";
    i(7):="one,  the  step sizes are adjusted according to the difficulty of";
    i(8):="the current path being followed.  In this case, a difficult  path";
    i(9):="does not slow down the continuation for the other ones.          ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Number_of_Paths_Info;

  procedure Max_Steps_Info is

    i : array(1..3) of string(1..65);

  begin
    i(1):="  This parameter bounds the amount of work along a solution path.";
    i(2):="Continuation for the path being followed stops when the number of";
    i(3):="steps exceeds this threshold value.                              ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Max_Steps_Info;

  procedure Distance_Info is

    i : array(1..6) of string(1..65);

  begin
    i(1):="  Tracking solution paths is organized in two stages.  As t<1, no";
    i(2):="singularities  are  likely to occur and paths can be tracked more";
    i(3):="loosely.  When t gets close to 1, paths may converge to  singular";
    i(4):="solutions or diverge to infinity, so they have to be tracked more";
    i(5):="carefully.  This distance marks the  switch  towards  the  second";
    i(6):="stage, that is the so-called end game.                           ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Distance_Info;

  procedure Extrapolation_Order_Info is

    i : array(1..5) of string(1..65);

  begin
    i(1):="  The direction of a diverging solution path provides a normal to";
    i(2):="a  face  of the system that certificates the deficiency.  In this";
    i(3):="polyhedral end game, extrapolation is used to  determine  winding";
    i(4):="numbers.   When  the  order of the extrapolator equals zero, then";
    i(5):="the polyhedral end game is turned off.                           ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Extrapolation_Order_Info;

  procedure Re_Run_Info is

    i : array(1..4) of string(1..65);

  begin
    i(1):="  Solution paths may be clustered at the end of the continuation.";
    i(2):="If  the  number  of  re-runs  is  higher than one, then clustered";
    i(3):="solution paths  are  re-computed  by  following  the  paths  more";
    i(4):="closely.  Setting this value to zero turns off this facility.    ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Re_Run_Info;

  procedure Predictor_Info is

    i : array(1..11) of string(1..65);

  begin
    i( 1):="  We distinguish between predictors for the solution vector x and";
    i( 2):="the  continuation  parameter  t.   The  secant predictor for x is";
    i( 3):="based on linear extrapolation, whereas the tangent predictor uses";
    i( 4):="the Jacobian matrix to extrapolate.                              ";
    i( 5):="  For t, three  different  predictors  are  provided.   The  real";
    i( 6):="predictor  simply  augments  t  with  the step size.  The complex";
    i( 7):="predictor takes values in the complex plane to circumvent  values";
    i( 8):="of  t  for  which  the corresponding system is badly conditioned.";
    i( 9):="The geometric predictor  generates  values  for  t  so  that  the";
    i(10):="distances  to  the  target  value behave like a geometric series.";
    i(11):="The latter predictor is useful in a polyhedral end game.         ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Predictor_Info;

  procedure Min_Step_Size_Info is

    i : array(1..4) of string(1..65);

  begin
    i(1):="  This is a lower bound on the step  size  for  the  continuation";
    i(2):="parameter.  Continuation for the path being followed stops when a";
    i(3):="required reduction of the step size yields a step length that  is";
    i(4):="lower than the minimum step size.                                ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Min_Step_Size_Info;

  procedure Max_Step_Size_Info is

    i : array(1..4) of string(1..65);

  begin
    i(1):="  This is an upper bound on the step size  for  the  continuation";
    i(2):="parameter.   An  expansion  of  the  step size that yields a step";
    i(3):="length that is larger than this threshold sets the step  size  to";
    i(4):="this threshold.                                                  ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Max_Step_Size_Info;

  procedure Reduction_Factor_Info is

    i : array(1..3) of string(1..65);

  begin
    i(1):="  When the corrector  does  not  manage  to  reach  the  required";
    i(2):="accuracy  within  the allowed number of iterations, the step size";
    i(3):="will be reduced by multiplication with the reduction factor.     ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Reduction_Factor_Info;

  procedure Expansion_Factor_Info is

    i : array(1..4) of string(1..65);

  begin
    i(1):="  When the corrector has reached the required accuracy within the";
    i(2):="allowed  number  of iterations, for a number of consecutive times";
    i(3):="higher than the  expansion  threshold,  the  step  size  will  be";
    i(4):="enlarged by multiplication by the expansion factor.              ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Expansion_Factor_Info;

  procedure Expansion_Threshold_Info is

    i : array(1..7) of string(1..65);

  begin
    i(1):="  Before enlarging the step size along  a  path,  the  number  of";
    i(2):="consecutive  successful  steps must pass the expansion threshold.";
    i(3):="A step is successful when the  corrector  manages  to  reach  the";
    i(4):="required accuracy within the allowed number of iterations.       ";
    i(5):="  At the end of the path, this threshold  bounds  the  number  of";
    i(6):="consecutive  equal  guesses  of the winding number $m$ before the";
    i(7):="polyhedral end game can modify it.                               ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Expansion_Threshold_Info;

  procedure Number_of_Iterations_Info is

    i : array(1..7) of string(1..65);

  begin
    i(1):="  The corrector stops when the desired  accuracy  is  reached  or";
    i(2):="when it has exhausted its maximum number of iterations.          ";
    i(3):="  A low maximum enforces quadratic convergence and keeps the path";
    i(4):="tracker close to the solution paths.                             ";
    i(5):="  A higher maximum may be needed at the end of the solution path,";
    i(6):="when  quadratic  convergence  can  no  longer  be obtained due to";
    i(7):="singularities.                                                   ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Number_of_Iterations_Info;

  procedure Relative_Residual_Info is

    i : array(1..5) of string(1..65);

  begin
    i(1):="  The  residual  is  the  norm  of  the  vector  obtained   after";
    i(2):="evaluating  the  current  approximation  vector in the polynomial";
    i(3):="system.  The corrector stops when the  residual  divided  by  the";
    i(4):="norm  of  the  approximation vector is lower than or equal to the";
    i(5):="required precision.                                              ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Relative_Residual_Info;

  procedure Absolute_Residual_Info is

    i : array(1..4) of string(1..65);

  begin
    i(1):="  The  residual  is  the  norm  of  the  vector  obtained   after";
    i(2):="evaluating  the  current  approximation  vector in the polynomial";
    i(3):="system.  The corrector stops when the residual is lower  than  or";
    i(4):="equal to the required precision.                                 ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Absolute_Residual_Info;

  procedure Relative_Correction_Info is

    i : array(1..4) of string(1..65);

  begin
    i(1):="  The correction is the norm of the last vector  used  to  update";
    i(2):="the  current  approximation vector.  The corrector stops when the";
    i(3):="correction divided by the norm of  the  approximation  vector  is";
    i(4):="lower than or equal to the required precision.                   ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Relative_Correction_Info;

  procedure Absolute_Correction_Info is

    i : array(1..3) of string(1..65);

  begin
    i(1):="  The correction is the norm of the last vector  used  to  update";
    i(2):="the  current  approximation vector.  The corrector stops when the";
    i(3):="correction is lower than or equal to the required precision.     ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Absolute_Correction_Info;

  procedure Singular_Info is

    i : array(1..2) of string(1..65);

  begin
    i(1):="  A condition is considered as singular when the condition number";
    i(2):="of the Jacobian matrix passes the given threshold value.         ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Singular_Info;

  procedure Clustered_Info is

    i : array(1..3) of string(1..65);

  begin
    i(1):="  Two solutions are considered as  clustered  when  the  distance";
    i(2):="between  all  corresponding  components is smaller than the given";
    i(3):="threshold value.                                                 ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Clustered_Info;

  procedure Infinity_Info is

    i : array(1..6) of string(1..65);

  begin
    i(1):="  A solution is considered to diverge to infinity when  its  norm";
    i(2):="passes  the given threshold value, in case of affine coordinates,";
    i(3):="or, in case of projective coordinates, when the added  coordinate";
    i(4):="becomes   smaller  than  the  inverse  of  the  threshold  value.";
    i(5):="Continuation for the path being followed stops when  it  diverges";
    i(6):="to infinity.                                                     ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Infinity_Info;

  procedure Exit_Info is

    i : constant string
      :="By typing 0, the current continuation parameters are used.";

  begin
    put_line(i);
  end Exit_Info;

  procedure Menu_Predictor ( pred_type : out natural32 ) is

  -- DESCRIPTION :
  --   Shows the menu with available predictors and lets the user choose.

    ans : character;
    m : array(0..10) of string(1..65);

  begin
    m(0):="MENU with Predictors for Solution x & Continuation Parameter t : ";
    m(1):="  0. Secant  : linear extrapolation & Real      increment for t; ";
    m(2):="  1. Secant  :  based on 2 previous & Complex   increment for t; ";
    m(3):="  2. Secant  :    solution values   & Geometric series    for t; ";
    m(4):="  3. Tangent : linear extrapolation & Real      increment for t; ";
    m(5):="  4. Tangent :   based on Jacobi    & Complex   increment for t; ";
    m(6):="  5. Tangent :  matrix at solution  & Geometric series    for t; ";
    m(7):="  6. Hermite : cubic extrapolation  & Real      increment for t; ";
    m(8):="  7. Quadratic : quadratic extrapol & Real      increment for t; ";
    m(9):="  8. Cubic   : cubic extrapolation  & Real      increment for t; ";
   m(10):="  9. Quartic :   quartic extrapol   & Real      increment for t; ";
    for i in m'range loop
      put_line(m(i));
    end loop;
    put("Select predictor by typing number between 0 and 9 : ");
    Ask_Alternative(ans,"0123456789");
    pred_type := Characters_and_Numbers.Convert(ans);
  end Menu_Predictor;

-- TARGET ROUTINES :

  procedure put is
  begin
    Continuation_Parameters_io.put(Standard_Output);
  end put;

  procedure put ( file : in file_type ) is

    banner_path : constant string := Predictor_Banner(predictor_path_type);
    banner_endg : constant string := Predictor_Banner(predictor_endg_type);

  begin

    put_line(file,"GLOBAL MONITOR : ");
    put(file,"  1. the condition of the homotopy           : ");
    put(file,condition,1); new_line(file);
    put(file,"  2. number of paths tracked simultaneously  : ");
    put(file,block_size,1); new_line(file);
    put(file,"  3. maximum number of steps along a path    : "); 
    put(file,max_steps,1); new_line(file);
    put(file,"  4. distance from target to start end game  : "); 
    put(file,start_end_game,1,3,3); new_line(file);
    put(file,"  5. order of extrapolator in end game       : ");
    put(file,endext_order,1); new_line(file);
    put(file,"  6. maximum number of re-runs               : ");
    put(file,max_reruns,1); new_line(file);

    put_line(file,"STEP CONTROL (PREDICTOR) : "
                          & "                   along path : end game");

    put(file,"  7: 8. type " & banner_path & ":"
                                         & banner_endg & " : ");
    put(file,predictor_path_type,1);    put(file,"         : ");
    put(file,predictor_endg_type,1); new_line(file);
    put(file,"  9:10. minimum step size                    : ");
    put(file,min_path_step_size,1,3,3);         put(file," : ");
    put(file,min_endg_step_size,1,3,3); new_line(file);
    put(file," 11:12. maximum step size                    : ");
    put(file,max_path_step_size,1,3,3);         put(file," : ");
    put(file,max_endg_step_size,1,3,3); new_line(file);
    put(file," 13:14. reduction factor for step size       : ");
    put(file,reduction_path_factor,1,3,3);      put(file," : ");
    put(file,reduction_endg_factor,1,3,3); new_line(file);
    put(file," 15:16. expansion factor for step size       : ");
    put(file,expansion_path_factor,1,3,3);      put(file," : ");
    put(file,expansion_endg_factor,1,3,3); new_line(file);
    put(file," 17:18. expansion threshold                  : ");
    put(file,success_path_steps,1);     put(file,"         : ");
    put(file,success_endg_steps,1); new_line(file);

    put_line(file,"PATH CLOSENESS (CORRECTOR) : " 
                            & "                 along path : end game");

    put(file," 19:20. maximum number of iterations         : ");
    put(file,max_path_iter,1);          put(file,"         : ");
    put(file,max_endg_iter,1); new_line(file);
    put(file," 21:22. relative precision for residuals     : ");
    put(file,relative_path_residual,1,3,3);     put(file," : ");
    put(file,relative_endg_residual,1,3,3); new_line(file);
    put(file," 23:24. absolute precision for residuals     : ");
    put(file,absolute_path_residual,1,3,3);     put(file," : ");
    put(file,absolute_endg_residual,1,3,3); new_line(file);
    put(file," 25:26. relative precision for corrections   : ");
    put(file,relative_path_correction,1,3,3);   put(file," : ");
    put(file,relative_endg_correction,1,3,3); new_line(file);
    put(file," 27:28. absolute precision for corrections   : ");
    put(file,absolute_path_correction,1,3,3);   put(file," : ");
    put(file,absolute_endg_correction,1,3,3); new_line(file);

    put_line(file,"SOLUTION TOLERANCES : "
                     & "                        along path : end game");

    put(file," 29:30. inverse condition of Jacobian        : ");
    put(file,tol_path_inverse_condition,1,3,3); put(file," : ");
    put(file,tol_endg_inverse_condition,1,3,3); new_line(file);
    put(file," 31:32. clustering of solutions              : ");
    put(file,tol_path_distance,1,3,3);          put(file," : ");
    put(file,tol_endg_distance,1,3,3); new_line(file);
    put(file," 33:34. solution at infinity                 : ");
    put(file,tol_path_at_infinity,1,3,3);       put(file," : ");
    put(file,tol_endg_at_infinity,1,3,3); new_line(file);

  end put;

  procedure get_parameter ( k : in natural32 ) is
  begin
    case k is
-- GLOBAL MONITOR :
      when  1 => put("  condition parameter : ");  Read_Natural(condition);
                 Continuation_Parameters.Tune(condition);
      when  2 => put("  number of paths tracked simultaneously : ");
                 Read_Positive(integer(block_size));
      when  3 => put("  maximum number of steps along a path : ");
                 Read_Positive(integer(max_steps));
      when  4 => put("  distance from target to start end game : ");
                 Read_Positive_Float(start_end_game);
      when  5 => put("  order of extrapolator in end game : ");
                 Read_Natural(endext_order);
                 Continuation_Parameters.Tune_Endgm_Pred(endext_order);
      when  6 => put("  maximum number of re-runs : ");
                 Read_Natural(max_reruns);
-- STEP CONTROL (PREDICTOR) :
      when  7 => Menu_Predictor(predictor_path_type);
      when  8 => Menu_Predictor(predictor_endg_type);
      when  9 => put("  minimum step size along the path : ");
                 Read_Positive_Float(min_path_step_size);
      when 10 => put("  minimum step size at end of path : ");
                 Read_Positive_Float(min_endg_step_size);
      when 11 => put("  maximum step size along the path : ");
                 Read_Positive_Float(max_path_step_size);
      when 12 => put("  maximum step size at end of path : ");
                 Read_Positive_Float(max_endg_step_size);
      when 13 => put("  reduction factor for step size along the path : ");
                 Read_Positive_Float(reduction_path_factor);
      when 14 => put("  reduction factor for step size at end of path : ");
                 Read_Positive_Float(reduction_endg_factor);
      when 15 => put("  expansion factor for step size along the path : ");
                 Read_Positive_Float(expansion_path_factor);
      when 16 => put("  expansion factor for step size at end of path : ");
                 Read_Positive_Float(expansion_endg_factor);
      when 17 => put("  expansion threshold along the path : ");
                 Read_Positive(integer(success_path_steps));
      when 18 => put("  expansion threshold at end of path : ");
                 Read_Positive(integer(success_endg_steps));
-- PATH CLOSENESS (CORRECTOR) :
      when 19 => put("  maximum number of iterations along the path : ");
                 Read_Positive(integer(max_path_iter));
      when 20 => put("  maximum number of iterations at end of path : ");
                 Read_Positive(integer(max_endg_iter));
      when 21 => put("  relative precision for residual along the path : ");
                 Read_Positive_Float(relative_path_residual);
      when 22 => put("  relative precision for residual at end of path : ");
                 Read_Positive_Float(relative_endg_residual);
      when 23 => put("  absolute precision for residual along the path : ");
                 Read_Positive_Float(absolute_path_residual);
      when 24 => put("  absolute precision for residual at end of path : ");
                 Read_Positive_Float(absolute_endg_residual);
      when 25 => put("  relative precision for correction along the path : ");
                 Read_Positive_Float(relative_path_correction);
      when 26 => put("  relative precision for correction at end of path : ");
                 Read_Positive_Float(relative_endg_correction);
      when 27 => put("  absolute precision for correction along the path : ");
                 Read_Positive_Float(absolute_path_correction);
      when 28 => put("  absolute precision for correction at end of path : ");
                 Read_Positive_Float(absolute_endg_correction);
-- SOLUTION TOLERANCES :
      when 29 => put("  inverse condition of Jacobian along the path : ");
                 Read_Positive_Float(tol_path_inverse_condition);
      when 30 => put("  inverse condition of Jacobian at end of path : ");
                 Read_Positive_Float(tol_endg_inverse_condition);
      when 31 => put("  clustering of solutions along the path : ");
                 Read_Positive_Float(tol_path_distance);
      when 32 => put("  clustering of solutions at end of path : ");
                 Read_Positive_Float(tol_endg_distance);
      when 33 => put("  solution at infinity along the path : ");
                 Read_Positive_Float(tol_path_at_infinity);
      when 34 => put("  solution at infinity at end of path : ");
                 Read_Positive_Float(tol_endg_at_infinity);
      when others => null;
    end case;
  end get_parameter;

  procedure get_parameter_with_info ( k : in natural32 ) is

    ans : character;

  begin
    new_line;
    case k is
      when  1 => Condition_Info;                             -- GLOBAL MONITOR 
      when  2 => Number_of_Paths_Info;
      when  3 => Max_Steps_Info;
      when  4 => Distance_Info;
      when  5 => Extrapolation_Order_Info;
      when  6 => Re_Run_Info;
      when  7 |  8 => Predictor_Info;              -- STEP CONTROL (PREDICTOR)
      when  9 | 10 => Min_Step_Size_Info;
      when 11 | 12 => Max_Step_Size_Info;
      when 13 | 14 => Reduction_Factor_Info;
      when 15 | 16 => Expansion_Factor_Info;
      when 17 | 18 => Expansion_Threshold_Info;
      when 19 | 20 => Number_of_Iterations_Info; -- PATH CLOSENESS (CORRECTOR)
      when 21 | 22 => Relative_Residual_Info;
      when 23 | 24 => Absolute_Residual_Info;
      when 25 | 26 => Relative_correction_Info;
      when 27 | 28 => Absolute_Correction_Info;
      when 29 | 30 => Singular_Info;                    -- SOLUTION TOLERANCES 
      when 31 | 32 => Clustered_Info;
      when 33 | 34 => Infinity_Info;
      when others  => Exit_Info;
    end case;
    new_line;
    if k /= 0 then
      put("Do you want to change the value of this parameter ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        get_parameter(k);
      end if;
    end if;
  end get_parameter_with_info;

  procedure get ( k : out natural32 ) is

    tmp : string(1..10);
    cnt : natural;
    ans : character;
    nb : integer32 := -1;
    valid,info,done : boolean;

  begin
    loop
      put("Type a number to change (0 to exit), preceded by i for info : ");
      get_line(tmp,cnt);
      valid := true; info := false;
      for i in 1..cnt loop
        if tmp(i) = 'i' then
          info := true;
        elsif tmp(i) /= ' '
             and Characters_and_Numbers.Convert(tmp(i)) = 10 then
          valid := false;
        end if;
        exit when not valid;
      end loop;
      if not valid then
        put_line("Invalid answer.  Please try again.");
      else
        done := false;
        loop
          declare
            nb_n : natural32;
          begin
            nb_n := Characters_and_Numbers.Convert(tmp(1..cnt));
            nb := integer32(nb_n);
            done := true;
          exception
            when others => 
              if nb < 0 then
                put_line("No natural number.  Please try again.");
                put("The current string is :"); put_line(tmp(1..cnt));
              end if;
          end;
          exit when done;
        end loop;
      end if;
      exit when (nb >= 0);
    end loop;
    if info then
      get_parameter_with_info(natural32(nb));
      if nb = 0 then
        put("Do you want to leave this menu ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'n'
         then nb := 1;
        end if;
      end if;
    else
      get_parameter(natural32(nb));
    end if;
    k := natural32(nb);
  end get;

end Continuation_Parameters_io;
