with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Test_Standard_Coeff_Homotopy;
with Test_DoblDobl_Coeff_Homotopy;

procedure ts_evalhomt is

-- DESCRIPTION :
--   Development of the evaluation of (1-t)*f + t*g
--   where monomials of f and g are shared
--   as in a coefficient-parameter homotopy.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and runs tests.

    ans : character;

  begin
    new_line;
    put_line("MENU for the precision : ");
    put_line("  1. double precision");
    put_line("  2. double double precision");
    put("Type 1 or 2 to select the precision : ");
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Test_Standard_Coeff_Homotopy.Main;
      when '2' => Test_DoblDobl_Coeff_Homotopy.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_evalhomt;
