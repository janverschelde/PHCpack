with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
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

end Pack_Continuation_Parameters;
