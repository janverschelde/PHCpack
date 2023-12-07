with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Standard_Speel_Convolutions;
with Test_DoblDobl_Speel_Convolutions;
with Test_TripDobl_Speel_Convolutions;
with Test_QuadDobl_Speel_Convolutions;
with Test_PentDobl_Speel_Convolutions;
with Test_OctoDobl_Speel_Convolutions;
with Test_DecaDobl_Speel_Convolutions;
with Test_HexaDobl_Speel_Convolutions;

procedure ts_speelser is

-- DESCRIPTION :
--   A vectorized version of the computation of Speelpenning products
--   for truncated power series.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree, the dimension,
  --   and the precision of the coefficients.

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
      when '0' => Test_Standard_Speel_Convolutions.Main;
      when '1' => Test_DoblDobl_Speel_Convolutions.Main;
      when '2' => Test_TripDobl_Speel_Convolutions.Main;
      when '3' => Test_QuadDobl_Speel_Convolutions.Main;
      when '4' => Test_PentDobl_Speel_Convolutions.Main;
      when '5' => Test_OctoDobl_Speel_Convolutions.Main;
      when '6' => Test_DecaDobl_Speel_Convolutions.Main;
      when '7' => Test_HexaDobl_Speel_Convolutions.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_speelser;
