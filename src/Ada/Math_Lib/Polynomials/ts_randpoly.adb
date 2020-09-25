with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Test_Standard_Random_Systems;
with Test_DoblDobl_Random_Systems;
with Test_TripDobl_Random_Systems;
with Test_QuadDobl_Random_Systems;
with Test_PentDobl_Random_Systems;
with Test_OctoDobl_Random_Systems;
with Test_DecaDobl_Random_Systems;
with Test_Multprec_Random_Systems;

procedure ts_randpoly is

-- DESCRIPTION :
--   Test on random polynomials.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and runs tests.

    ans : character;

  begin
    new_line;
    put_line("Generation of random dense and sparse polynomial systems.");
    new_line;
    put_line("MENU for the precision : ");
    put_line("  1. double precision");
    put_line("  2. double double precision");
    put_line("  3. triple double precision");
    put_line("  4. quad double precision");
    put_line("  5. penta double precision");
    put_line("  6. octo double precision");
    put_line("  7. deca double precision");
    put_line("  8. arbitrary multiprecision");
    put("Type 1, 2, 3, 4, 5, 6, 7, or 8 to make a choice : ");
    Ask_Alternative(ans,"12345678");
    case ans is
      when '1' => Test_Standard_Random_Systems.Main;
      when '2' => Test_DoblDobl_Random_Systems.Main;
      when '3' => Test_TripDobl_Random_Systems.Main;
      when '4' => Test_QuadDobl_Random_Systems.Main;
      when '5' => Test_PentDobl_Random_Systems.Main;
      when '6' => Test_OctoDobl_Random_Systems.Main;
      when '7' => Test_DecaDobl_Random_Systems.Main;
      when '8' => Test_Multprec_Random_Systems.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_randpoly;
