with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_SD_Newton_Matrix_Series;
with Test_DD_Newton_Matrix_Series;
with Test_TD_Newton_Matrix_Series;
with Test_QD_Newton_Matrix_Series;
with Test_PD_Newton_Matrix_Series;
with Test_OD_Newton_Matrix_Series;
with Test_DA_Newton_Matrix_Series;
with Test_HD_Newton_Matrix_Series;

procedure ts_sernew is

-- DESCRIPTION :
--   Tests the application of Newton's method to compute series solutions.

  precision : character;

begin
  new_line;
  put_line("MENU for the working precision :");
  put_line("  0. double precision");
  put_line("  1. double double precision");
  put_line("  2. triple double precision");
  put_line("  3. quad double precision");
  put_line("  4. penta double precision");
  put_line("  5. octo double precision");
  put_line("  6. deca double precision");
  put_line("  7. hexa double precision");
  put("Type 0, 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
  Ask_Alternative(precision,"01234567");
  case precision is 
    when '0' => Test_SD_Newton_Matrix_Series.Main;
    when '1' => Test_DD_Newton_Matrix_Series.Main;
    when '2' => Test_TD_Newton_Matrix_Series.Main;
    when '3' => Test_QD_Newton_Matrix_Series.Main;
    when '4' => Test_PD_Newton_Matrix_Series.Main;
    when '5' => Test_OD_Newton_Matrix_Series.Main;
    when '6' => Test_DA_Newton_Matrix_Series.Main;
    when '7' => Test_HD_Newton_Matrix_Series.Main;
    when others => null;
  end case;
end ts_sernew;
