with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Numbers_io;

package body QuadDobl_Quad_Parameters is

-- RESET to DEFAULTS :

  procedure Reset is
  begin
    max_step_size := create(0.1);
    reduction_multiplier := create(0.5);
    expansion_multiplier := create(1.5);
    expansion_threshold := 3;
    increment_tolerance := create(1.0E-8);
    residual_tolerance := create(1.0E-8);
    max_corrector_steps := 4;
    max_predictor_steps := 500;
    determinant_tolerance := create(1.0E-12);
  end Reset;

-- REPORTING and TUNING :

  procedure Show is
  begin
    Show(standard_output);
  end Show;

  procedure Show ( file : in file_type ) is
  begin
    put_line(file,"***** CURRENT SWEEP PARAMETER SETTINGS *****");
    put(file,"  1. maximal number of predictor steps : ");
    put(file,max_predictor_steps,1); new_line(file);
    put(file,"  2. maximal step size for predictor   : ");
    put(file,max_step_size,3); new_line(file);
    put(file,"  3. multiplier to reduce step size    : ");
    put(file,reduction_multiplier,3); new_line(file);
    put(file,"  4. multiplier to expand step size    : ");
    put(file,expansion_multiplier,3); new_line(file);
    put(file,"  5. threshold to delay expansion      : ");
    put(file,expansion_threshold,1); new_line(file);
    put(file,"  6. tolerance on evaluation residual  : ");
    put(file,residual_tolerance,3); new_line(file);
    put(file,"  7. tolerance on smallest increment   : "); 
    put(file,increment_tolerance,3); new_line(file);
    put(file,"  8. maximal number of corrector steps : ");
    put(file,max_corrector_steps,1); new_line(file);
    put(file,"  9. tolerance on Jacobian determinant : ");
    put(file,determinant_tolerance,3); new_line(file);
  end Show;

  procedure Tune is

    ans : character;
 
  begin
    loop
      Show;
      put("Type a number to change settings, R to reset, or 0 to exit : ");
      Ask_Alternative(ans,"0123456789R");
      exit when (ans = '0');
      case ans is
        when '1' =>
          put("Give new maximal number of predictor steps : ");
          Numbers_io.Read_Positive(integer(max_predictor_steps));
        when '2' =>
          put("Give new maximal step size for predictor : ");
          Numbers_io.Read_Positive_Quad_Double(max_step_size);
        when '3' =>
          put("Give new multiplier to reduce the step size : ");
          Numbers_io.Read_Positive_Quad_Double(reduction_multiplier);
        when '4' => 
          put("Give new multiplier to expand the step size : ");
          Numbers_io.Read_Positive_Quad_Double(expansion_multiplier);
        when '5' =>
          put("Give new threshold to delay the expansion : ");
          Numbers_io.Read_Positive(integer(expansion_threshold));
        when '6' =>
          put("Give new tolerance of evaluation residual : ");
          Numbers_io.Read_Positive_Quad_Double(residual_tolerance);
        when '7' =>
          put("Give new tolerance on smallest increment : ");
          Numbers_io.Read_Positive_Quad_Double(increment_tolerance);
        when '8' =>
          put("Give new maximal number of corrector steps : ");
          Numbers_io.Read_Positive(integer(max_corrector_steps));
        when '9' =>
          put("Give new tolerance of Jacobian determinant : ");
          Numbers_io.Read_Positive_Quad_Double(determinant_tolerance);
        when others => Reset;
      end case;
    end loop;
  end Tune;

end QuadDobl_Quad_Parameters;
