with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Cells_Interface;

function use_celcon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Cells_Interface;

  begin
    case job is
      when 0 => return Cells_Read_Floating_Mixed_Cells(vrblvl);
      when 1 => return Cells_Write_Floating_Mixed_Cells(vrblvl);
      when 2 => return Cells_Number_of_Floating_Mixed_Cells(a,vrblvl);
      when 3 => return Cells_Dimension_of_Floating_Mixed_Cells(a,vrblvl);
      when 4 => return Cells_Get_Floating_Mixture(a,b,vrblvl);
      when 5 => return Cells_Floating_Supports_Size(a,b,vrblvl);
      when 6 => return Cells_Get_Floating_Support_Point(a,b,c,vrblvl);
      when 7 => return Cells_Floating_Normal(a,c,vrblvl);
      when 8 => return Cells_Floating_Cell_Size(a,b,vrblvl);
      when 9 => return Cells_Get_Floating_Cell_Point(a,b,c,vrblvl);
      when 10 => return Cells_Floating_Mixed_Volume(a,b,vrblvl);
      when 11 => return Cells_Set_Floating_Mixture(a,b,vrblvl);
      when 12 => return Cells_Add_Floating_Support_Point(a,b,c,vrblvl);
      when 13 => return Cells_Add_Floating_Mixed_Cell(a,b,c,vrblvl);
      when 14 => return Cells_Floating_Clear(vrblvl);
      when 15 => return Cells_Get_Floating_Mixed_Cell(a,b,c,vrblvl);
      when 16 => return Cells_Make_Standard_Coefficient_System(vrblvl);
      when 17 => return Cells_Read_Standard_Coefficient_System(vrblvl);
      when 18 => return Cells_Write_Standard_Coefficient_System(vrblvl);
      when 19 => return Cells_Standard_System_into_Container(vrblvl);
      when 20 => return Cells_Standard_System_from_Container(vrblvl);
      when 21 => return Cells_Standard_Polyhedral_Homotopy(vrblvl);
      when 22 => return Cells_Standard_Start_Solve(a,b,vrblvl);
      when 23 => return Cells_Standard_Track_One_Path(a,b,vrblvl);
      when 24 => return Cells_Standard_TarSol_into_Container(a,b,vrblvl);
      when 25 => return Cells_Standard_Permute(vrblvl);
      when 26 => return Cells_Make_DoblDobl_Coefficient_System(vrblvl);
      when 27 => return Cells_Read_DoblDobl_Coefficient_System(vrblvl);
      when 28 => return Cells_Write_DoblDobl_Coefficient_System(vrblvl);
      when 29 => return Cells_DoblDobl_System_into_Container(vrblvl);
      when 30 => return Cells_DoblDobl_System_from_Container(vrblvl);
      when 31 => return Cells_DoblDobl_Polyhedral_Homotopy(vrblvl);
      when 32 => return Cells_DoblDobl_Start_Solve(a,b,vrblvl);
      when 33 => return Cells_DoblDobl_Track_One_Path(a,b,vrblvl);
      when 34 => return Cells_DoblDobl_TarSol_into_Container(a,b,vrblvl);
      when 35 => return Cells_DoblDobl_Permute(vrblvl);
      when 36 => return Cells_Make_QuadDobl_Coefficient_System(vrblvl);
      when 37 => return Cells_Read_QuadDobl_Coefficient_System(vrblvl);
      when 38 => return Cells_Write_QuadDobl_Coefficient_System(vrblvl);
      when 39 => return Cells_QuadDobl_System_into_Container(vrblvl);
      when 40 => return Cells_QuadDobl_System_from_Container(vrblvl);
      when 41 => return Cells_QuadDobl_Polyhedral_Homotopy(vrblvl);
      when 42 => return Cells_QuadDobl_Start_Solve(a,b,vrblvl);
      when 43 => return Cells_QuadDobl_Track_One_Path(a,b,vrblvl);
      when 44 => return Cells_QuadDobl_TarSol_into_Container(a,b,vrblvl);
      when 45 => return Cells_QuadDobl_Permute(vrblvl);
      when 46 => return Cells_Floating_Mixed_Volume(a,vrblvl);
      when 47 => return Cells_Set_Floating_Number_of_Supports(a,vrblvl);
      when 48 => return Cells_Standard_StaSol_into_Container(a,b,vrblvl);
      when 49 => return Cells_DoblDobl_StaSol_into_Container(a,b,vrblvl);
      when 50 => return Cells_QuadDobl_StaSol_into_Container(a,b,vrblvl);
      when 51 => return Cells_Read_Integer_Mixed_Cells(vrblvl);
      when 52 => return Cells_Write_Integer_Mixed_Cells(vrblvl);
      when 53 => return Cells_Number_of_Integer_Mixed_Cells(a,vrblvl);
      when 54 => return Cells_Dimension_of_Integer_Mixed_Cells(a,vrblvl);
      when 55 => return Cells_Get_Integer_Mixture(a,b,vrblvl);
      when 56 => return Cells_Integer_Supports_Size(a,b,vrblvl);
      when 57 => return Cells_Get_Integer_Support_Point(a,b,c,vrblvl);
      when 58 => return Cells_Integer_Normal(a,c,vrblvl);
      when 59 => return Cells_Integer_Cell_Size(a,b,vrblvl);
      when 60 => return Cells_Get_Integer_Cell_Point(a,b,c,vrblvl);
      when 61 => return Cells_Integer_Mixed_Volume(a,b,vrblvl);
      when 62 => return Cells_Set_Integer_Number_of_Supports(a,vrblvl);
      when 63 => return Cells_Set_Integer_Mixture(a,b,vrblvl);
      when 64 => return Cells_Add_Integer_Support_Point(a,b,c,vrblvl);
      when 65 => return Cells_Add_Integer_Mixed_Cell(a,b,c,vrblvl);
      when 66 => return Cells_Integer_Clear(vrblvl);
      when 67 => return Cells_Get_Integer_Mixed_Cell(a,b,c,vrblvl);
      when 68 => return Cells_Make_Integer_Subdivision(vrblvl);
      when 69 => return Cells_Is_Stable(a,vrblvl);
      when 70 => return Cells_Number_of_Original_Cells(a,vrblvl);
      when 71 => return Cells_Number_of_Stable_Cells(a,vrblvl);
      when 72 => return Cells_Standard_Stable_Solve(a,b,vrblvl);
      when 73 => return Cells_DoblDobl_Stable_Solve(a,b,vrblvl);
      when 74 => return Cells_QuadDobl_Stable_Solve(a,b,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  exception
    when others => put("Exception raised in use_celcon handling job ");
                   put(job,1); put_line(".  Will try to ignore...");
                   return 1;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_celcon;
