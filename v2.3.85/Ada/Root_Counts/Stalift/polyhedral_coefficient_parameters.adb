with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Numbers_io;                         use Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

package body Polyhedral_Coefficient_Parameters is

  procedure Write ( file : in file_type ) is
  begin
    put_line(file,
     "****************** CURRENT CONTINUATION PARAMETERS *****************");
    put_line(file,"GLOBAL :");
    put(file,"  1. maximal #steps along a path      : ");
                                    put(file,max_nb_steps,2); new_line(file);
    put_line(file,"STEP CONTROL (PREDICTOR) :                             ");
    put(file,"  2. start for continuation parameter :");
                                    put(file,min_infinity,3); new_line(file);
    put(file,"  3. maximum step size for predictor  :");
                                   put(file,max_pred_step,3); new_line(file);
    put(file,"  4. reduction factor for step size   :");
                                put(file,reduction_factor,3); new_line(file);
    put(file,"  5. expansion factor for step size   :");
                                put(file,expansion_factor,3); new_line(file);
    put(file,"  6. expansion threshold              : ");
                             put(file,expansion_threshold,2); new_line(file);
    put_line(file,"PATH CLOSENESS (CORRECTOR) :                           ");
    put(file,"  7. maximum #corrector iterations    : ");
                                   put(file,max_corr_iter,2); new_line(file);
    put(file,"  8. tolerance on root accuracy       :");
                                        put(file,tol_root,3); new_line(file);
    put_line(file,
     "********************************************************************");
  end Write;

  procedure Tune is

    ans : character;

  begin
    loop
      Write(standard_output);
      put("Type a number to change a parameter (0 to exit) : ");
      Ask_Alternative(ans,"012345678");
      exit when ans = '0';
      case ans is
        when '1' => put("Give new maximum #steps along path : ");
                    Read_Natural(max_nb_steps);
        when '2' => put("Give new start for continuation parameter : ");
                    Read_Double_Float(min_infinity);
        when '3' => put("Give new maximum step size : ");
                    Read_Double_Float(max_pred_step);
        when '4' => put("Give new reduction factor : ");
                    Read_Double_Float(reduction_factor);
        when '5' => put("Give new expansion factor : ");
                    Read_Double_Float(expansion_factor);
        when '6' => put("Give new expansion threshold : ");
                    Read_Natural(expansion_threshold);
        when '7' => put("Give new #corrector iterations : ");
                    Read_Natural(max_corr_iter);
        when '8' => put("Give new tolerance on root accuracy : ");
                    Read_Double_Float(tol_root);
        when others => null;
      end case;
    end loop;
  end Tune;

end Polyhedral_Coefficient_Parameters;
