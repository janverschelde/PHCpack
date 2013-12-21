with integer_io;                         use integer_io;
with File_Scanning;                      use File_Scanning;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Continuation_Parameters;            use Continuation_Parameters;
with Process_io;                         use Process_io;

package body Scanners_for_Continuation is

-- AUXILIAIRIES :

  procedure Scan_Natural ( file : in file_type; n : in out natural ) is

  -- DESCRIPTION :
  --   Scans the current line of the file for a natural number,
  --   which should occur after a `:'.

    found : boolean;

  begin
    Scan(file,":",found); get(file,n); skip_line(file);
  exception
    when others => put_line("INCORRECTLY SPECIFIED NATURAL NUMBER");
                   raise;
  end Scan_Natural;

  procedure Scan_Float ( file : in file_type; f : in out double_float ) is

  -- DESCRIPTION :
  --   Scans the current line of the file for a floating point number,
  --   which should occur after a `:'.

    found : boolean;

  begin
    Scan(file,":",found); get(file,f); skip_line(file);
  exception
    when others => put_line("INCORRECTLY SPECIFIED FLOATING POINT NUMBER");
                   raise;
  end Scan_Float;

  procedure Scan_Complex ( file : in file_type; c : in out Complex_Number ) is

  -- DESCRIPTION :
  --   Scans the current line of the file for a complex number,
  --   which should occur after a `:'.

    found : boolean;

  begin
    Scan(file,":",found); get(file,c); skip_line(file);
  exception
    when others => put_line("INCORRECTLY SPECIFIED COMPLEX NUMBER");
                   raise;
  end Scan_Complex;

-- TARGET PROCEDURES :

  procedure Scan_Homotopy_Parameters
              ( file : in file_type;
                k : out natural; a : out Complex_Number ) is

    found : boolean;
    kk : natural := 0;
    aa : Complex_Number := Create(0.0);

  begin
    Scan_and_Skip(file,"HOMOTOPY PARAMETERS",found);
    if found
     then Scan_Natural(file,kk); Scan_Complex(file,aa);
    end if;
    k := kk; a := aa;
  exception
    when others => put_line("INCORRECTLY SPECIFIED HOMOTOPY PARAMETERS");
                   raise;
  end Scan_Homotopy_Parameters;

  procedure Scan_Continuation_Parameters ( file : in file_type ) is

    found : boolean;

  begin
    Scan_and_Skip(file,"CONTINUATION PARAMETERS",found);
    if found
     then null;

    --   Scan_and_Skip(file,"CONDITION",found);
    --   Scan_Natural(file,condition);

    --   Scan_and_Skip(file,"MONITOR",found);
    --   Scan_Natural(file,block_size);
    --   Scan_Natural(file,predictor_type);
    --   Scan_Natural(file,max_steps);
    --   Scan_Float(file,start_end_game);
    --   Scan_Natural(file,max_reruns);

    --   Scan_and_Skip(file,"STEP CONTROL",found);
    --   Scan_Float(file,min_step_size);
    --   Scan_Float(file,max_step_size);
    --   Scan_Float(file,reduction_factor);
    --   Scan_Float(file,expansion_factor);
    --   Scan_Natural(file,success_steps);

    --   Scan_and_Skip(file,"PATH CLOSENESS",found);
    --   Scan_Natural(file,max_iter_path);
    --   Scan_Natural(file,max_iter_end);
    --   Scan_Float(file,relative_path_residual);
    --   Scan_Float(file,absolute_path_residual);
    --   Scan_Float(file,relative_path_correction);
    --   Scan_Float(file,absolute_path_correction);
    --   Scan_Float(file,relative_end_residual);
    --   Scan_Float(file,absolute_end_residual);
    --   Scan_Float(file,relative_end_correction);
    --   Scan_Float(file,absolute_end_correction);

    --   Scan_and_Skip(file,"SOLUTIONS",found);
    --   Scan_Float(file,tol_inverse_condition);
    --   Scan_Float(file,tol_distance);
    --   Scan_Float(file,tol_at_infinity);

    end if;
  exception
    when others => put_line("INCORRECTLY SPECIFIED CONTINUATION PARAMETERS");
                   raise;
  end Scan_Continuation_Parameters;

  procedure Scan_Output_Parameter ( file : in file_type; op : out natural ) is

    found : boolean;
    oc : natural := 0;

  begin
    Scan(file,"OUTPUT PARAMETER",found); 
    if found
     then Scan_Natural(file,oc);
          case oc is
            when 1 => Set_output_code(s);
            when 2 => Set_output_code(p);
            when 3 => Set_output_code(c);
            when 4 => Set_output_code(sp);
            when 5 => Set_output_code(sc);
            when 6 => Set_output_code(pc);
            when 7 => Set_output_code(spc);
            when others => Set_output_code(nil);
          end case;
     else Set_output_code(nil);
    end if;
    op := oc;
  exception
    when others => put_line("INCORRECTLY SPECIFIED OUTPUT PARAMETER");
                   raise;
  end Scan_Output_Parameter;

end Scanners_for_Continuation;
