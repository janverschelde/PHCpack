with text_io;                           use text_io;
with File_Management_Interface;
with Standard_PolySys_Interface;
with DoblDobl_PolySys_Interface;
with QuadDobl_PolySys_Interface;
with Linear_Products_Interface;
with Newton_Interface;
with Path_Trackers_Interface;
with Cascade_Homotopy_Interface;
with Diagonal_Homotopy_Interface;

function use_track ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer;
                     vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use File_Management_Interface;
    use Standard_PolySys_Interface;
    use DoblDobl_PolySys_Interface;
    use QuadDobl_PolySys_Interface;
    use Linear_Products_Interface;
    use Newton_Interface;
    use Path_Trackers_Interface;
    use Cascade_Homotopy_Interface;
    use Diagonal_Homotopy_Interface;

  begin
    case job is
      when -1 => return Newton_Standard_Polynomial_Refine(b,c,vrblvl-1);
      when 0 => return Standard_PolySys_Prompt_for_Target(vrblvl-1);
      when 1 => return Standard_PolySys_Prompt_for_Start(vrblvl-1);
      when 2 => return Path_Trackers_Standard_Homotopy_Random(vrblvl-1);
      when 3 => return Path_Trackers_Standard_Homotopy_Gamma(a,c,vrblvl-1);
      when 4 => return Path_Trackers_Standard_Homotopy_Clear(vrblvl-1);
      when 5 => return Path_Trackers_Standard_Silent_Track(a,b,c,vrblvl-1);
      when 6 => return Path_Trackers_Standard_Report_Track(a,b,c,vrblvl-1);
      when 7 => return Path_Trackers_Standard_Write_Solution(a,b,c,vrblvl-1);
      when 8 => return File_Management_Write_String(a,b,vrblvl-1);
      when 9 => return File_Management_Write_Integers(a,b,vrblvl-1);
      when 10 => return File_Management_Write_Doubles(a,c,vrblvl-1);
      when 11 => return Standard_PolySys_Read_Target_on_File(a,b,vrblvl-1);
      when 12 => return Standard_PolySys_Read_Start_on_File(a,b,vrblvl-1);
      when 13 => return Linear_Products_System_Read(a,b,vrblvl-1);
      when 14 => return Cascade_Homotopy_Standard_Polynomial(vrblvl-1);
      when 15 =>
        return Diagonal_Homotopy_Standard_Polynomial_Make(a,b,vrblvl-1);
      when 16 => return Diagonal_Homotopy_Prompt_Set(a,b,vrblvl-1);
      when 17 => return Diagonal_Homotopy_Reset_Input(a,b,vrblvl-1);
      when 18 => return Diagonal_Homotopy_Cascade_Dimension(a,b,vrblvl-1);
      when 19 => return Diagonal_Homotopy_Standard_Hyperset(a,b,vrblvl-1);
      when 20 => return Diagonal_Homotopy_Standard_Collapse(a,vrblvl-1);
      when 21 => return Cascade_Homotopy_Cut_Slack(a,vrblvl-1);
     -- tracking in double double precision :
      when 22 => return Path_Trackers_DoblDobl_Homotopy_Random(vrblvl-1);
      when 23 => return Path_Trackers_DoblDobl_Homotopy_Gamma(a,c,vrblvl-1);
      when 24 => return Path_Trackers_DoblDobl_Homotopy_Clear(vrblvl-1);
      when 25 => return Path_Trackers_DoblDobl_Silent_Track(a,b,c,vrblvl-1);
      when 26 => return Path_Trackers_DoblDobl_Report_Track(a,b,c,vrblvl-1);
      when 27 => return Path_Trackers_DoblDobl_Write_Solution(a,b,c,vrblvl-1);
      when 28 => return Cascade_Homotopy_DoblDobl_Polynomial(vrblvl-1);
     -- tracking in quad double precision :
      when 32 => return Path_Trackers_QuadDobl_Homotopy_Random(vrblvl-1);
      when 33 => return Path_Trackers_QuadDobl_Homotopy_Gamma(a,c,vrblvl-1);
      when 34 => return Path_Trackers_QuadDobl_Homotopy_Clear(vrblvl-1);
      when 35 => return Path_Trackers_QuadDobl_Silent_Track(a,b,c,vrblvl-1);
      when 36 => return Path_Trackers_QuadDobl_Report_Track(a,b,c,vrblvl-1);
      when 37 => return Path_Trackers_QuadDobl_Write_Solution(a,b,c,vrblvl-1);
      when 38 => return Cascade_Homotopy_QuadDobl_Polynomial(vrblvl-1);
     -- redefining diagonal homotopies ...
      when 40 =>
        return Diagonal_Homotopy_Standard_Polynomial_Set(a,b,vrblvl-1);
      when 41 =>
        return Diagonal_Homotopy_Standard_Start_Solutions(a,b,vrblvl-1);
      when 42 => return Diagonal_Homotopy_Symbols_Doubler(a,b,vrblvl-1);
     -- diagonal homotopy in double double and quad double precision
      when 43 =>
        return Diagonal_Homotopy_DoblDobl_Polynomial_Make(a,b,vrblvl-1);
      when 44 =>
        return Diagonal_Homotopy_QuadDobl_Polynomial_Make(a,b,vrblvl-1);
      when 45 =>
        return Diagonal_Homotopy_DoblDobl_Start_Solutions(a,b,vrblvl-1);
      when 46 =>
        return Diagonal_Homotopy_QuadDobl_Start_Solutions(a,b,vrblvl-1);
      when 47 => return Diagonal_Homotopy_DoblDobl_Collapse(a,vrblvl-1);
      when 48 => return Diagonal_Homotopy_QuadDobl_Collapse(a,vrblvl-1);
     -- double double and quad double witness sets for hypersurface
      when 49 =>
        return Diagonal_Homotopy_DoblDobl_Polynomial_Set(a,b,vrblvl-1);
      when 50 =>
        return Diagonal_Homotopy_QuadDobl_Polynomial_Set(a,b,vrblvl-1);
     -- multiprecision versions to create homotopy :
      when 52 => return Path_Trackers_Multprec_Homotopy_Random(vrblvl-1);
      when 53 => return Path_Trackers_Multprec_Homotopy_Gamma(a,c,vrblvl-1);
      when 54 => return Path_Trackers_Multprec_Homotopy_Clear(vrblvl-1);
     -- crude path trackers
      when 55 => return Path_Trackers_Standard_Crude_Track(a,vrblvl-1);
      when 56 => return Path_Trackers_DoblDobl_Crude_Track(a,vrblvl-1);
      when 57 => return Path_Trackers_QuadDobl_Crude_Track(a,vrblvl-1);
     -- homotopies for Laurent systems
      when 58 => return Cascade_Homotopy_Standard_Laurent(vrblvl-1);
      when 59 => return Cascade_Homotopy_DoblDobl_Laurent(vrblvl-1);
      when 60 => return Cascade_Homotopy_QuadDobl_Laurent(vrblvl-1);
      when 61 => return Diagonal_Homotopy_Standard_Laurent_Make(a,b,vrblvl-1);
      when 62 => return Diagonal_Homotopy_DoblDobl_Laurent_Make(a,b,vrblvl-1);
      when 63 => return Diagonal_Homotopy_QuadDobl_Laurent_Make(a,b,vrblvl-1);
      when 64 =>
        return Diagonal_Homotopy_Standard_Laurential_Set(a,b,vrblvl-1);
      when 65 =>
        return Diagonal_Homotopy_DoblDobl_Laurential_Set(a,b,vrblvl-1);
      when 66 =>
        return Diagonal_Homotopy_QuadDobl_Laurential_Set(a,b,vrblvl-1);
      when 67 => return DoblDobl_PolySys_Prompt_for_Target(vrblvl-1);
      when 68 => return QuadDobl_PolySys_Prompt_for_Target(vrblvl-1);
     -- more homotopy operations on Laurent systems
      when 69 =>
        return Path_Trackers_Standard_Laurent_Homotopy_Gamma(a,c,vrblvl-1);
      when 70 =>
        return Path_Trackers_DoblDobl_Laurent_Homotopy_Gamma(a,c,vrblvl-1);
      when 71 =>
        return Path_Trackers_QuadDobl_Laurent_Homotopy_Gamma(a,c,vrblvl-1);
      when 72 =>
        return Path_Trackers_Standard_Laurent_Homotopy_Clear(vrblvl-1);
      when 73 =>
        return Path_Trackers_DoblDobl_Laurent_Homotopy_Clear(vrblvl-1);
      when 74 =>
        return Path_Trackers_QuadDobl_Laurent_Homotopy_Clear(vrblvl-1);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_track;
