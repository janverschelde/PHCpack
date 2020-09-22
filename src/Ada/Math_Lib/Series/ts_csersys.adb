with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Test_Standard_CSeries_Systems;
with Test_DoblDobl_CSeries_Systems;
with Test_QuadDobl_CSeries_Systems;

procedure ts_csersys is

-- DESCRIPTION :
--   Tests the methods on systems of series polynomials.

  procedure Main is

  -- DESCRIPTION :
  --   Displays the menu for the working precision and
  --   calls the test corresponding to the choice made.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Test_Standard_CSeries_Systems.Main;
      when '1' => Test_DoblDobl_CSeries_Systems.Main;
      when '2' => Test_QuadDobl_CSeries_Systems.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_csersys;
