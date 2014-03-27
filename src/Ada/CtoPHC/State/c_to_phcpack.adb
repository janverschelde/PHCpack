with text_io;                          use text_io;
with PHCpack_Operations;               use PHCpack_Operations;
with PHCpack_Operations_io;

function C_to_PHCpack ( job : integer32;
                        number_of_tasks : natural32 ) return integer32 is

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
    put_line(" 28. clear the quad double data.");
  end Write_Menu;

begin
  case job is
    when 0 => Write_Menu;
   -- operations on standard doubles :
    when 1 => PHCpack_Operations_io.Read_Target_System;
    when 2 => PHCpack_Operations_io.Write_Target_System;
    when 3 => PHCpack_Operations_io.Read_Start_System;
    when 4 => PHCpack_Operations_io.Write_Start_System;
    when 5 => PHCpack_Operations_io.Write_Start_Solutions;
    when 6 => return Solve_by_Standard_Homotopy_Continuation(number_of_tasks);
    when 7 => PHCpack_Operations_io.Write_Target_Solutions;
    when 8 => Standard_Clear;
    when 9 => Define_Output_File;
   -- operations on double doubles :
    when 11 => PHCpack_Operations_io.Read_DoblDobl_Target_System;
    when 12 => PHCpack_Operations_io.Write_DoblDobl_Target_System;
    when 13 => PHCpack_Operations_io.Read_DoblDobl_Start_System;
    when 14 => PHCpack_Operations_io.Write_DoblDobl_Start_System;
    when 15 => PHCpack_Operations_io.Write_DoblDobl_Start_Solutions;
    when 16 => return Solve_by_DoblDobl_Homotopy_Continuation(number_of_tasks);
    when 17 => PHCpack_Operations_io.Write_DoblDobl_Target_Solutions;
    when 18 => DoblDobl_Clear;
   -- operations on quad doubles :
    when 21 => PHCpack_Operations_io.Read_QuadDobl_Target_System;
    when 22 => PHCpack_Operations_io.Write_QuadDobl_Target_System;
    when 23 => PHCpack_Operations_io.Read_QuadDobl_Start_System;
    when 24 => PHCpack_Operations_io.Write_QuadDobl_Start_System;
    when 25 => PHCpack_Operations_io.Write_QuadDobl_Start_Solutions;
    when 26 => return Solve_by_QuadDobl_Homotopy_Continuation(number_of_tasks);
    when 27 => PHCpack_Operations_io.Write_QuadDobl_Target_Solutions;
    when 28 => QuadDobl_Clear;
    when others => put_line("  Sorry, this operation is not defined.");
                   return 1;
  end case;
  return 0;
end C_to_PHCpack;
