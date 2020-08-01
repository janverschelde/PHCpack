with text_io;                           use text_io;
with Step_Trackers_Interface;

function use_nxtsol ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Step_Trackers_Interface;

  begin
    case job is
      when 0 => return Step_Trackers_Standard_Homotopy(a,c,vrblvl);
      when 1 => return Step_Trackers_DoblDobl_Homotopy(a,c,vrblvl);
      when 2 => return Step_Trackers_QuadDobl_Homotopy(a,c,vrblvl);
      when 3 => return Step_Trackers_Set_Standard_Solution(a,vrblvl);
      when 4 => return Step_Trackers_Set_DoblDobl_Solution(a,vrblvl);
      when 5 => return Step_Trackers_Set_QuadDobl_Solution(a,vrblvl);
      when 6 => return Step_Trackers_Next_Standard_Solution(a,vrblvl);
      when 7 => return Step_Trackers_Next_DoblDobl_Solution(a,vrblvl);
      when 8 => return Step_Trackers_Next_QuadDobl_Solution(a,vrblvl);
      when 9 => return Step_Trackers_Standard_Clear(vrblvl);
      when 10 => return Step_Trackers_DoblDobl_Clear(vrblvl);
      when 11 => return Step_Trackers_QuadDobl_Clear(vrblvl);
      when 12 => return Step_Trackers_Multprec_Homotopy(a,b,vrblvl);
      when 13 => return Step_Trackers_Set_Multprec_Solution(a,vrblvl);
      when 14 => return Step_Trackers_Next_Multprec_Solution(a,vrblvl);
      when 15 => return Step_Trackers_Multprec_Clear(vrblvl);
      when 16 => return Step_Trackers_Varbprec_Homotopy(a,b,vrblvl);
      when 17 => return Step_Trackers_Set_Varbprec_Solution(a,b,vrblvl);
      when 18 => return Step_Trackers_Next_Varbprec_Solution(a,b,vrblvl);
      when 19 => return Step_Trackers_Varbprec_Clear(vrblvl);
      when 20 => return Step_Trackers_Get_Varbprec_Solution(a,b,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_nxtsol;
