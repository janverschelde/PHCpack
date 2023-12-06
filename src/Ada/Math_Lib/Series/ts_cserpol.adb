with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Test_Standard_CSeries_Polynomials;
with Test_DoblDobl_CSeries_Polynomials;
with Test_TripDobl_CSeries_Polynomials;
with Test_QuadDobl_CSeries_Polynomials;
with Test_PentDobl_CSeries_Polynomials;
with Test_OctoDobl_CSeries_Polynomials;
with Test_DecaDobl_CSeries_Polynomials;
with Test_HexaDobl_CSeries_Polynomials;

procedure ts_cserpol is

-- DESCRIPTION :
--   Tests the development of polynomials in several variables,
--   with truncated power series as coefficients.
 
  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and then calls
  --   the main test in that precision.

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
    put("Type 0, 1, 2, 3, 4, 5, 6, or 7 to select the working precision : ");
    Ask_Alternative(prc,"01234567");
    case prc is 
      when '0' => Test_Standard_CSeries_Polynomials.Main;
      when '1' => Test_DoblDobl_CSeries_Polynomials.Main;
      when '2' => Test_TripDobl_CSeries_Polynomials.Main;
      when '3' => Test_QuadDobl_CSeries_Polynomials.Main;
      when '4' => Test_PentDobl_CSeries_Polynomials.Main;
      when '5' => Test_OctoDobl_CSeries_Polynomials.Main;
      when '6' => Test_DecaDobl_CSeries_Polynomials.Main;
      when '7' => Test_HexaDobl_CSeries_Polynomials.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cserpol;
