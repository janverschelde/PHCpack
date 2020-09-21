with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Test_Standard_CSeries_Polynomials;
with Test_DoblDobl_CSeries_Polynomials;
with Test_TripDobl_CSeries_Polynomials;
with Test_QuadDobl_CSeries_Polynomials;
with Test_PentDobl_CSeries_Polynomials;
with Test_OctoDobl_CSeries_Polynomials;
with Test_DecaDobl_CSeries_Polynomials;

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
    put_line("  1. double precision");
    put_line("  2. double double precision");
    put_line("  3. triple double precision");
    put_line("  4. quad double precision");
    put_line("  5. penta double precision");
    put_line("  6. octo double precision");
    put_line("  7. deca double precision");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the working precision : ");
    Ask_Alternative(prc,"1234567");
    case prc is 
      when '1' => Test_Standard_CSeries_Polynomials.Main;
      when '2' => Test_DoblDobl_CSeries_Polynomials.Main;
      when '3' => Test_TripDobl_CSeries_Polynomials.Main;
      when '4' => Test_QuadDobl_CSeries_Polynomials.Main;
      when '5' => Test_PentDobl_CSeries_Polynomials.Main;
      when '6' => Test_OctoDobl_CSeries_Polynomials.Main;
      when '7' => Test_DecaDobl_CSeries_Polynomials.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cserpol;
