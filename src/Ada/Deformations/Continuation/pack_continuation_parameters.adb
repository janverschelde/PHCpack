with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Continuation_Parameters;

package body Pack_Continuation_Parameters is

  function Get return Standard_Floating_Vectors.Vector is

    v : Standard_Floating_Vectors.Vector(1..34);

  begin
    v(1) := double_float(Continuation_Parameters.condition);
    v(2) := double_float(Continuation_Parameters.block_size);
    v(3) := double_float(Continuation_Parameters.max_steps);
    v(4) := Continuation_Parameters.start_end_game;
    v(5) := double_float(Continuation_Parameters.endext_order);
    v(6) := double_float(Continuation_Parameters.max_reruns);
    v(7) := double_float(Continuation_Parameters.predictor_path_type);
    v(8) := double_float(Continuation_Parameters.predictor_endg_type);
    v(9) := Continuation_Parameters.min_path_step_size;
    v(10) := Continuation_Parameters.min_endg_step_size;
    v(11) := Continuation_Parameters.max_path_step_size;
    v(12) := Continuation_Parameters.max_endg_step_size;
    v(13) := Continuation_Parameters.reduction_path_factor;
    v(14) := Continuation_Parameters.reduction_endg_factor;
    v(15) := Continuation_Parameters.expansion_path_factor;
    v(16) := Continuation_Parameters.expansion_endg_factor;
    v(17) := double_float(Continuation_Parameters.success_path_steps);
    v(18) := double_float(Continuation_Parameters.success_endg_steps);
    v(19) := double_float(Continuation_Parameters.max_path_iter);
    v(20) := double_float(Continuation_Parameters.max_endg_iter);
    v(21) := Continuation_Parameters.relative_path_residual;
    v(22) := Continuation_Parameters.relative_endg_residual;
    v(23) := Continuation_Parameters.absolute_path_residual;
    v(24) := Continuation_Parameters.absolute_endg_residual;
    v(25) := Continuation_Parameters.relative_path_correction;
    v(26) := Continuation_Parameters.relative_endg_correction;
    v(27) := Continuation_Parameters.absolute_path_correction; 
    v(28) := Continuation_Parameters.absolute_endg_correction;
    v(29) := Continuation_Parameters.tol_path_inverse_condition;
    v(30) := Continuation_Parameters.tol_endg_inverse_condition;
    v(31) := Continuation_Parameters.tol_path_distance;
    v(32) := Continuation_Parameters.tol_endg_distance;
    v(33) := Continuation_Parameters.tol_path_at_infinity;
    v(34) := Continuation_Parameters.tol_endg_at_infinity;
    return v;
  end Get;

  function Get_Value ( k : natural32 ) return double_float is

    res : double_float;

  begin
    case k is
      when 1 => res := double_float(Continuation_Parameters.condition);
      when 2 => res := double_float(Continuation_Parameters.block_size);
      when 3 => res := double_float(Continuation_Parameters.max_steps);
      when 4 => res := Continuation_Parameters.start_end_game;
      when 5 => res := double_float(Continuation_Parameters.endext_order);
      when 6 => res := double_float(Continuation_Parameters.max_reruns);
      when 7 =>
        res := double_float(Continuation_Parameters.predictor_path_type);
      when 8 =>
        res := double_float(Continuation_Parameters.predictor_endg_type);
      when 9 => res := Continuation_Parameters.min_path_step_size;
      when 10 => res := Continuation_Parameters.min_endg_step_size;
      when 11 => res := Continuation_Parameters.max_path_step_size;
      when 12 => res := Continuation_Parameters.max_endg_step_size;
      when 13 => res := Continuation_Parameters.reduction_path_factor;
      when 14 => res := Continuation_Parameters.reduction_endg_factor;
      when 15 => res := Continuation_Parameters.expansion_path_factor;
      when 16 => res := Continuation_Parameters.expansion_endg_factor;
      when 17 =>
        res := double_float(Continuation_Parameters.success_path_steps);
      when 18 =>
        res := double_float(Continuation_Parameters.success_endg_steps);
      when 19 => res := double_float(Continuation_Parameters.max_path_iter);
      when 20 => res := double_float(Continuation_Parameters.max_endg_iter);
      when 21 => res := Continuation_Parameters.relative_path_residual;
      when 22 => res := Continuation_Parameters.relative_endg_residual;
      when 23 => res := Continuation_Parameters.absolute_path_residual;
      when 24 => res := Continuation_Parameters.absolute_endg_residual;
      when 25 => res := Continuation_Parameters.relative_path_correction;
      when 26 => res := Continuation_Parameters.relative_endg_correction;
      when 27 => res := Continuation_Parameters.absolute_path_correction; 
      when 28 => res := Continuation_Parameters.absolute_endg_correction;
      when 29 => res := Continuation_Parameters.tol_path_inverse_condition;
      when 30 => res := Continuation_Parameters.tol_endg_inverse_condition;
      when 31 => res := Continuation_Parameters.tol_path_distance;
      when 32 => res := Continuation_Parameters.tol_endg_distance;
      when 33 => res := Continuation_Parameters.tol_path_at_infinity;
      when 34 => res := Continuation_Parameters.tol_endg_at_infinity;
      when others => res := -1.0; 
    end case;
    return res;
  end Get_Value;

  procedure Write ( v : in Standard_Floating_Vectors.Vector ) is
  begin
    Write(standard_output,v);
  end Write;

  procedure Write ( file : in file_type;
                    v : in Standard_Floating_Vectors.Vector ) is
  begin
    put(file,natural32(v(1)),1); put(file," ");
    put(file,natural32(v(2)),1); put(file," ");
    put(file,natural32(v(3)),1); put(file," ");
    put(file,v(4),2); put(file," ");
    put(file,natural32(v(5)),1); put(file," ");
    put(file,natural32(v(6)),1); put(file," ");
    put(file,natural32(v(7)),1); put(file," ");
    put(file,natural32(v(8)),1); put(file," "); new_line(file);
    for i in 9..16 loop
      put(file,v(integer32(i)),2);
    end loop;
    new_line(file);
    put(file,natural32(v(17)),1); put(file," ");
    put(file,natural32(v(18)),1); put(file," ");
    put(file,natural32(v(19)),1); put(file," ");
    put(file,natural32(v(20)),1); put(file," ");
    for i in 21..27 loop
      put(file,v(integer32(i)),2);
    end loop;
    new_line(file);
    for i in 28..34 loop
      put(file,v(integer32(i)),2);
    end loop;
    new_line(file);
  end Write;

  procedure Set ( v : in Standard_Floating_Vectors.Vector ) is
  begin
    Continuation_Parameters.condition := natural32(v(1));
    Continuation_Parameters.block_size := natural32(v(2));
    Continuation_Parameters.max_steps := natural32(v(3));
    Continuation_Parameters.start_end_game := v(4);
    Continuation_Parameters.endext_order := natural32(v(5));
    Continuation_Parameters.max_reruns := natural32(v(6));
    Continuation_Parameters.predictor_path_type := natural32(v(7));
    Continuation_Parameters.predictor_endg_type := natural32(v(8));
    Continuation_Parameters.min_path_step_size := v(9);
    Continuation_Parameters.min_endg_step_size := v(10);
    Continuation_Parameters.max_path_step_size := v(11);
    Continuation_Parameters.max_endg_step_size := v(12);
    Continuation_Parameters.reduction_path_factor := v(13);
    Continuation_Parameters.reduction_endg_factor := v(14);
    Continuation_Parameters.expansion_path_factor := v(15);
    Continuation_Parameters.expansion_endg_factor := v(16);
    Continuation_Parameters.success_path_steps := natural32(v(17));
    Continuation_Parameters.success_endg_steps := natural32(v(18));
    Continuation_Parameters.max_path_iter := natural32(v(19));
    Continuation_Parameters.max_endg_iter := natural32(v(20));
    Continuation_Parameters.relative_path_residual := v(21);
    Continuation_Parameters.relative_endg_residual := v(22);
    Continuation_Parameters.absolute_path_residual := v(23);
    Continuation_Parameters.absolute_endg_residual := v(24);
    Continuation_Parameters.relative_path_correction := v(25);
    Continuation_Parameters.relative_endg_correction := v(26);
    Continuation_Parameters.absolute_path_correction := v(27);
    Continuation_Parameters.absolute_endg_correction := v(28);
    Continuation_Parameters.tol_path_inverse_condition := v(29);
    Continuation_Parameters.tol_endg_inverse_condition := v(30);
    Continuation_Parameters.tol_path_distance := v(31);
    Continuation_Parameters.tol_endg_distance := v(32);
    Continuation_Parameters.tol_path_at_infinity := v(33);
    Continuation_Parameters.tol_endg_at_infinity := v(34);
  end Set;

  procedure Set_Value ( k : in natural32; v : in double_float ) is
  begin
    case k is
      when  1 => Continuation_Parameters.condition := natural32(v);
      when  2 => Continuation_Parameters.block_size := natural32(v);
      when  3 => Continuation_Parameters.max_steps := natural32(v);
      when  4 => Continuation_Parameters.start_end_game := v;
      when  5 => Continuation_Parameters.endext_order := natural32(v);
      when  6 => Continuation_Parameters.max_reruns := natural32(v);
      when  7 => Continuation_Parameters.predictor_path_type := natural32(v);
      when  8 => Continuation_Parameters.predictor_endg_type := natural32(v);
      when  9 => Continuation_Parameters.min_path_step_size := v;
      when 10 => Continuation_Parameters.min_endg_step_size := v;
      when 11 => Continuation_Parameters.max_path_step_size := v;
      when 12 => Continuation_Parameters.max_endg_step_size := v;
      when 13 => Continuation_Parameters.reduction_path_factor := v;
      when 14 => Continuation_Parameters.reduction_endg_factor := v;
      when 15 => Continuation_Parameters.expansion_path_factor := v;
      when 16 => Continuation_Parameters.expansion_endg_factor := v;
      when 17 => Continuation_Parameters.success_path_steps := natural32(v);
      when 18 => Continuation_Parameters.success_endg_steps := natural32(v);
      when 19 => Continuation_Parameters.max_path_iter := natural32(v);
      when 20 => Continuation_Parameters.max_endg_iter := natural32(v);
      when 21 => Continuation_Parameters.relative_path_residual := v;
      when 22 => Continuation_Parameters.relative_endg_residual := v;
      when 23 => Continuation_Parameters.absolute_path_residual := v;
      when 24 => Continuation_Parameters.absolute_endg_residual := v;
      when 25 => Continuation_Parameters.relative_path_correction := v;
      when 26 => Continuation_Parameters.relative_endg_correction := v;
      when 27 => Continuation_Parameters.absolute_path_correction := v;
      when 28 => Continuation_Parameters.absolute_endg_correction := v;
      when 29 => Continuation_Parameters.tol_path_inverse_condition := v;
      when 30 => Continuation_Parameters.tol_endg_inverse_condition := v;
      when 31 => Continuation_Parameters.tol_path_distance := v;
      when 32 => Continuation_Parameters.tol_endg_distance := v;
      when 33 => Continuation_Parameters.tol_path_at_infinity := v;
      when 34 => Continuation_Parameters.tol_endg_at_infinity := v;
      when others => null;
    end case;
  end Set_Value;

end Pack_Continuation_Parameters;
