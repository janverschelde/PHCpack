with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Standard_NewtConvSteps;
with Test_DoblDobl_NewtConvSteps;
with Test_TripDobl_NewtConvSteps;
with Test_QuadDobl_NewtConvSteps;
with Test_PentDobl_NewtConvSteps;
with Test_OctoDobl_NewtConvSteps;
with Test_DecaDobl_NewtConvSteps;
with Test_HexaDobl_NewtConvSteps;

procedure ts_sernewcnv is

-- DESCRIPTION :
--   Tests the linearized Newton's method for power series,
--   on convolution circuits.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and then runs tests.

    ans : character;

  begin
    new_line;
    put_line("Linearized Newton on power series with convolution circuits.");
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
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(ans,"1234567");
    case ans is
      when '0' => Test_Standard_NewtConvSteps.Main;
      when '1' => Test_DoblDobl_NewtConvSteps.Main;
      when '2' => Test_TripDobl_NewtConvSteps.Main;
      when '3' => Test_QuadDobl_NewtConvSteps.Main;
      when '4' => Test_PentDobl_NewtConvSteps.Main;
      when '5' => Test_OctoDobl_NewtConvSteps.Main;
      when '6' => Test_DecaDobl_NewtConvSteps.Main;
      when '7' => Test_HexaDobl_NewtConvSteps.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sernewcnv;
