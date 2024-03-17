with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Pade_Continuation_Interface;

function use_padcon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Pade_Continuation_Interface;

  begin
    case job is
      when 0 => return Pade_Continuation_Parameters_Set_Defaults(vrblvl);
      when 1 => return Pade_Continuation_Parameters_Clear(vrblvl);
      when 2 => return Pade_Continuation_Parameters_Get_Value(a,b,c,vrblvl);
      when 3 => return Pade_Continuation_Parameters_Set_Value(a,b,c,vrblvl);
      when 4 => return Pade_Continuation_Track_Paths(a,b,c,vrblvl);
      when 5 => return Pade_Continuation_Artificial_Homotopy(a,b,vrblvl);
      when 6 => return Pade_Continuation_Set_Start_Solution(a,b,vrblvl);
      when 7 => return Pade_Continuation_Next_Step(a,b,vrblvl);
      when 8 => return Pade_Continuation_Set_Solution(a,b,vrblvl);
      when 9 => return Pade_Continuation_Clear_Data(a,vrblvl);
      when 10 => return Pade_Continuation_Pole_Radius(a,c,vrblvl);
      when 11 => return Pade_Continuation_Closest_Pole(a,c,vrblvl);
      when 12 => return Pade_Continuation_T_Value(a,c,vrblvl);
      when 13 => return Pade_Continuation_Step_Size(a,c,vrblvl);
      when 14 => return Pade_Continuation_Series_Coefficient(a,b,c,vrblvl);
      when 15 => return Pade_Continuation_Pade_Coefficient(a,b,c,vrblvl);
      when 16 => return Pade_Continuation_Get_Pole(a,b,c,vrblvl);
      when 17 => return Pade_Continuation_Parameters_Write(vrblvl);
      when 18 => return Pade_Continuation_Natural_Homotopy(a,b,vrblvl);
      when 19 => return Pade_Continuation_Series_Step(a,c,vrblvl);
      when 20 => return Pade_Continuation_Pole_Step(a,c,vrblvl);
      when 21 => return Pade_Continuation_Estimated_Distance(a,c,vrblvl);
      when 22 => return Pade_Continuation_Hessian_Step(a,c,vrblvl);
      when 23 => return Pade_Continuation_Parameters_Reset_Values(a,vrblvl);
      when 24 => return Pade_Continuation_Set_Predicted_Solution(a,b,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return -1;
    end case;
  exception
    when others => put("Exception raised in use_padcon handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_padcon;
