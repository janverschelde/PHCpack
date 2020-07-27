with text_io;                           use text_io;
with Sweep_Interface;

function use_sweep ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer;
                     vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Sweep_Interface;

  begin
    case job is
      when 0 => return Sweep_Define_Parameters_Numerically(a,b,vrblvl);
      when 1 => return Sweep_Define_Parameters_Symbolically(a,b,vrblvl);
      when 2 => return Sweep_Number_of_Equations(a,vrblvl);
      when 3 => return Sweep_Number_of_Variables(a,vrblvl);
      when 4 => return Sweep_Number_of_Parameters(a,vrblvl);
      when 5 => return Sweep_Get_Parameters_Numerically(a,vrblvl);
      when 6 => return Sweep_Get_Parameters_Symbolically(a,b,vrblvl);
      when 7 => return Sweep_Parameters_Clear(vrblvl);
      when 8 => return Sweep_Set_Parameter_Values(a,b,c,vrblvl);
      when 9 => return Sweep_Get_Parameter_Values(a,c,vrblvl);
      when 10 => return Sweep_Complex_Convex_Parameter(a,c,vrblvl);
      when 11 => return Sweep_Real_Natural_Parameter(a,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_sweep;
