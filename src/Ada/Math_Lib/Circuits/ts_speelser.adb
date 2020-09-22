with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Test_Standard_Speel_Convolutions;
with Test_DoblDobl_Speel_Convolutions;
with Test_QuadDobl_Speel_Convolutions;

procedure ts_speelser is

-- DESCRIPTION :
--   A vectorized version of the computation of Speelpenning products
--   for truncated power series.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree, the dimension,
  --   and the precision of the coefficients.

    dim,deg : integer32 := 0;
    precision : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    case precision is
      when '0' => Test_Standard_Speel_Convolutions.Main;
      when '1' => Test_DoblDobl_Speel_Convolutions.Main;
      when '2' => Test_QuadDobl_Speel_Convolutions.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_speelser;
