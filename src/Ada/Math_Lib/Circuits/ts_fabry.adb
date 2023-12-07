with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Standard_Fabry;
with Test_DoblDobl_Fabry;
with Test_TripDobl_Fabry;
with Test_QuadDobl_Fabry;
with Test_PentDobl_Fabry;
with Test_OctoDobl_Fabry;
with Test_DecaDobl_Fabry;
with Test_HexaDobl_Fabry;

procedure ts_fabry is

-- DESCRIPTION :
--   Given a polynomial system, adds a parameter t to every coefficient
--   and runs the Newton's method on the power series.
--   The smallest ratio of the coefficients of the series will give
--   the convergence radius and the location for the nearest singularity.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the degree of the series,
  --   and for the precision.  Launches the tests.

    prc : character;

  begin
    new_line;
    put_line("Testing the Fabry ratio theorem ...");
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
    Ask_Alternative(prc,"01234567");
    case prc is
      when '0' => Test_Standard_Fabry.Main;
      when '1' => Test_DoblDobl_Fabry.Main;
      when '2' => Test_TripDobl_Fabry.Main;
      when '3' => Test_QuadDobl_Fabry.Main;
      when '4' => Test_PentDobl_Fabry.Main;
      when '5' => Test_OctoDobl_Fabry.Main;
      when '6' => Test_DecaDobl_Fabry.Main;
      when '7' => Test_HexaDobl_Fabry.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_fabry;
