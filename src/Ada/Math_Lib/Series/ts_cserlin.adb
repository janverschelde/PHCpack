with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Standard_Linearization;
with Test_DoblDobl_Linearization;
with Test_QuadDobl_Linearization;

procedure ts_cserlin is

-- DESCRIPTION :
--   Tests the linearization of solving linear systems of truncated series.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and then runs tests.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Test_Standard_Linearization.Main;
      when '1' => Test_DoblDobl_Linearization.Main;
      when '2' => Test_QuadDobl_Linearization.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cserlin;
