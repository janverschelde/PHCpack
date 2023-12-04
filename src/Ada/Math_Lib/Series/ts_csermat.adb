with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Test_Standard_Matrix_Series;
with Test_DoblDobl_Matrix_Series;
with Test_TripDobl_Matrix_Series;
with Test_QuadDobl_Matrix_Series;
with Test_PentDobl_Matrix_Series;
with Test_OctoDobl_Matrix_Series;
with Test_DecaDobl_Matrix_Series;
with Test_HexaDobl_Matrix_Series;

procedure ts_csermat is

-- DESCRIPTION :
--   Tests matrices of truncated dense power series.

  procedure Main is

  -- DESCRIPTION :
  --   Displays the menu of test procedures and runs the selected test.

    prc : character;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. triple double precision");
    put_line("  3. quad double precision");
    put_line("  4. penta double precision");
    put_line("  5. octo double precision");
    put_line("  6. deca double precision");
    put_line("  7. hexa double precision");
    put("Type 0, 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(prc,"01234567");
    case prc is
      when '0' => Test_Standard_Matrix_Series.Main;
      when '1' => Test_DoblDobl_Matrix_Series.Main;
      when '2' => Test_TripDobl_Matrix_Series.Main;
      when '3' => Test_QuadDobl_Matrix_Series.Main;
      when '4' => Test_PentDobl_Matrix_Series.Main;
      when '5' => Test_OctoDobl_Matrix_Series.Main;
      when '6' => Test_DecaDobl_Matrix_Series.Main;
      when '7' => Test_HexaDobl_Matrix_Series.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_csermat;
