with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Standard_Polynomials;
with Test_Standard_Laurentials;
with Test_Multprec_Polynomials;

procedure ts_poly is

-- DESCRIPTION :
--   Tests complex polynomials.

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for the precision and tests :");
    put_line("  0. test polynomials in standard double precision");
    put_line("  1. test Laurent polynomials in standard double precision");
    put_line("  2. test polynomials in arbitrary multiprecision");
    put("Type 0, 1, or 2 to select : "); Ask_Alternative(ans,"012");
    case ans is
      when '0' => Test_Standard_Polynomials.Main;
      when '1' => Test_Standard_Laurentials.Main;
      when '2' => Test_Multprec_Polynomials.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_poly;
