with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Test_Standard_CSeries_Polynomials;
with Test_DoblDobl_CSeries_Polynomials;
with Test_QuadDobl_CSeries_Polynomials;

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
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision;");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(prc,"012");
    case prc is 
      when '0' => Test_Standard_CSeries_Polynomials.Main;
      when '1' => Test_DoblDobl_CSeries_Polynomials.Main;
      when '2' => Test_QuadDobl_CSeries_Polynomials.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cserpol;
