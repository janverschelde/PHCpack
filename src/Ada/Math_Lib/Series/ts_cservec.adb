with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Standard_Vector_Series;
with Test_DoblDobl_Vector_Series;
with Test_TripDobl_Vector_Series;
with Test_QuadDobl_Vector_Series;
with Test_PentDobl_Vector_Series;
with Test_OctoDobl_Vector_Series;
with Test_DecaDobl_Vector_Series;

procedure ts_cservec is

-- DESCRIPTION :
--   Test on vectors of truncated power series.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for a precision and runs tests.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  1. standard double precision");
    put_line("  2. double double precision");
    put_line("  3. triple double precision");
    put_line("  4. quad double precision");
    put_line("  5. penta double precision");
    put_line("  6. octo double precision");
    put_line("  7. deca double precision");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(ans,"1234567");
    case ans is
      when '1' => Test_Standard_Vector_Series.Main;
      when '2' => Test_DoblDobl_Vector_Series.Main;
      when '3' => Test_TripDobl_Vector_Series.Main;
      when '4' => Test_QuadDobl_Vector_Series.Main;
      when '5' => Test_PentDobl_Vector_Series.Main;
      when '6' => Test_OctoDobl_Vector_Series.Main;
      when '7' => Test_DecaDobl_Vector_Series.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cservec;
