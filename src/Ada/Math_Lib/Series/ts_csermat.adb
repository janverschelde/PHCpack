with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Test_Standard_Matrix_Series;
with Test_DoblDobl_Matrix_Series;
with Test_QuadDobl_Matrix_Series;

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
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(prc,"012");
    case prc is
      when '0' => Test_Standard_Matrix_Series.Main;
      when '1' => Test_DoblDobl_Matrix_Series.Main;
      when '2' => Test_QuadDobl_Matrix_Series.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_csermat;
