with text_io;                           use text_io;
with Pieri_Interface;

function use_c2pieri ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Pieri_Interface;

  begin
    case job is
      when 0 => return Pieri_Write_Menu;
      when 1 => return Pieri_Initialize_Dimensions(a,vrblvl);
      when 2 => return Pieri_Initialize_Input_planes(a,b,c,vrblvl);
      when 3 => return Pieri_Initialize_Interpolation_points(a,b,c,vrblvl);
      when 4 => return Pieri_Store_Start_Pivots(a,b,vrblvl);
      when 5 => return Pieri_Store_Target_Pivots(a,b,vrblvl);
      when 6 => return Pieri_Store_Start_Coefficients(a,c,vrblvl);
      when 7 => return Pieri_Get_Target_Solution(a,c,vrblvl);
      when 8 => return Pieri_Silent_Track(vrblvl);
      when 9 => return Pieri_Report_Track(vrblvl);
      when 10 => return Pieri_Silent_Verify(c,vrblvl);
      when 11 => return Pieri_Report_Verify(c,vrblvl);
      when 12 => return Pieri_Clear(vrblvl);
      when 13 => return Pieri_Root_Count(a,b,vrblvl);
      when 14 => return Pieri_Localization_String(a,b,vrblvl);
      when 15 => return Pieri_Run_Homotopies(a,b,c,vrblvl);
      when 16 => return Pieri_Real_Osculating_Planes(a,c,vrblvl);
      when 17 => return Pieri_Make_Target_System(a,c,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2pieri;
