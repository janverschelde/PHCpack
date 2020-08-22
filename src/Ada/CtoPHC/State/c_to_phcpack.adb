with text_io;                          use text_io;
with PHCpack_Operations;               use PHCpack_Operations;
with PHCpack_Operations_io;

function C_to_PHCpack
           ( job : integer32;
             number_of_tasks : natural32;
             vrblvl : integer32 := 0 ) return integer32 is

  procedure Write_Menu is
  begin
    new_line;
    put_line("MENU to the operations in PHCpack :");
    put_line("  0. displays this menu;");
    put_line("  1. read target system;");
    put_line("  2. write target system;");
    put_line("  3. read start system;");
    put_line("  4. write start system;");
    put_line("  5. write start solutions;");
    put_line("  6. solve by standard homotopy continuation;");
    put_line("  7. write the target solutions;");
    put_line("  8. clear all data;");
    put_line("  9. define the output file.");
    put_line(" 11. read double double target system;");
    put_line(" 12. write double double target system;");
    put_line(" 13. read double double start system;");
    put_line(" 14. write double double start system;");
    put_line(" 15. write double double start solutions;");
    put_line(" 16. solve by double double homotopy continuation;");
    put_line(" 17. write the double double target solutions;");
    put_line(" 18. clear the double double data;");
    put_line(" 21. read quad double target system;");
    put_line(" 22. write quad double target system;");
    put_line(" 23. read quad double start system;");
    put_line(" 24. write quad double start system;");
    put_line(" 25. write quad double start solutions;");
    put_line(" 26. solve by quad double homotopy continuation;");
    put_line(" 27. write the quad double target solutions;");
    put_line(" 28. clear the quad double data;");
    put_line(" 29. read standard double start system;");
    put_line(" 30. write standard double start system;");
    put_line(" 31. read standard double target system;");
    put_line(" 32. write standard double target system;");
    put_line(" 33. read double double start system;");
    put_line(" 34. write double double start system;");
    put_line(" 35. read double double target system;");
    put_line(" 36. write double double target system;");
    put_line(" 37. read quad double start system;");
    put_line(" 38. write quad double start system;");
    put_line(" 39. read quad double target system;");
    put_line(" 40. write quad double target system;");
    put_line(" 41. clear the data for Laurent homotopy in double precision;");
    put_line(" 42. clear the data for Laurent homotopy with double doubles;");
    put_line(" 43. clear the data for Laurent homotopy with quad doubles.");
  end Write_Menu;

  function Handle_Jobs return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in C_to_PHCpack.Handle_Jobs ...");
    end if;
    case job is
      when 0 => Write_Menu;
     -- operations on standard doubles :
      when 1 => PHCpack_Operations_io.Read_Target_System;
      when 2 => PHCpack_Operations_io.Write_Target_System;
      when 3 => PHCpack_Operations_io.Read_Start_System;
      when 4 => PHCpack_Operations_io.Write_Start_System;
      when 5 => PHCpack_Operations_io.Write_Start_Solutions;
      when 6 =>
        return Solve_by_Standard_Homotopy_Continuation(number_of_tasks);
      when 7 => PHCpack_Operations_io.Write_Target_Solutions;
      when 8 => PHCpack_Operations.Standard_Clear;
      when 9 => Define_Output_File;
     -- operations on double doubles :
      when 11 => PHCpack_Operations_io.Read_DoblDobl_Target_System;
      when 12 => PHCpack_Operations_io.Write_DoblDobl_Target_System;
      when 13 => PHCpack_Operations_io.Read_DoblDobl_Start_System;
      when 14 => PHCpack_Operations_io.Write_DoblDobl_Start_System;
      when 15 => PHCpack_Operations_io.Write_DoblDobl_Start_Solutions;
      when 16 =>
        return Solve_by_DoblDobl_Homotopy_Continuation(number_of_tasks);
      when 17 => PHCpack_Operations_io.Write_DoblDobl_Target_Solutions;
      when 18 => PHCpack_Operations.DoblDobl_Clear;
     -- operations on quad doubles :
      when 21 => PHCpack_Operations_io.Read_QuadDobl_Target_System;
      when 22 => PHCpack_Operations_io.Write_QuadDobl_Target_System;
      when 23 => PHCpack_Operations_io.Read_QuadDobl_Start_System;
      when 24 => PHCpack_Operations_io.Write_QuadDobl_Start_System;
      when 25 => PHCpack_Operations_io.Write_QuadDobl_Start_Solutions;
      when 26 =>
        return Solve_by_QuadDobl_Homotopy_Continuation(number_of_tasks);
      when 27 => PHCpack_Operations_io.Write_QuadDobl_Target_Solutions;
      when 28 => PHCpack_Operations.QuadDobl_Clear;
      when 29 => PHCpack_Operations_io.Read_Start_Laurent_System;
      when 30 => PHCpack_Operations_io.Write_Start_Laurent_System;
      when 31 => PHCpack_Operations_io.Read_Target_Laurent_System;
      when 32 => PHCpack_Operations_io.Write_Target_Laurent_System;
      when 33 => PHCpack_Operations_io.Read_DoblDobl_Start_Laurent_System;
      when 34 => PHCpack_Operations_io.Write_DoblDobl_Start_Laurent_System;
      when 35 => PHCpack_Operations_io.Read_DoblDobl_Target_Laurent_System;
      when 36 => PHCpack_Operations_io.Write_DoblDobl_Target_Laurent_System;
      when 37 => PHCpack_Operations_io.Read_QuadDobl_Start_Laurent_System;
      when 38 => PHCpack_Operations_io.Write_QuadDobl_Start_Laurent_System;
      when 39 => PHCpack_Operations_io.Read_QuadDobl_Target_Laurent_System;
      when 40 => PHCpack_Operations_io.Write_QuadDobl_Target_Laurent_System;
      when 41 => PHCpack_Operations.Standard_Laurent_Clear;
      when 42 => PHCpack_Operations.DoblDobl_Laurent_Clear;
      when 43 => PHCpack_Operations.QuadDobl_Laurent_Clear;
      when 44 =>
        return Solve_by_Standard_Laurent_Homotopy_Continuation(number_of_tasks);
      when 45 =>
        return Solve_by_DoblDobl_Laurent_Homotopy_Continuation(number_of_tasks);
      when 46 =>
        return Solve_by_QuadDobl_Laurent_Homotopy_Continuation(number_of_tasks);
      when others => put_line("  Sorry, this operation is not defined.");
                     return 1;
    end case;
    return 0;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end C_to_PHCpack;
