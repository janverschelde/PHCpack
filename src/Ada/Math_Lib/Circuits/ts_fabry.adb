with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Standard_Fabry;
with Test_DoblDobl_Fabry;
with Test_TripDobl_Fabry;
with Test_QuadDobl_Fabry;
with Test_PentDobl_Fabry;
with Test_OctoDobl_Fabry;
with Test_DecaDobl_Fabry;

procedure ts_fabry is

-- DESCRIPTION :
--   Given a polynomial system, adds a parameter t to every coefficient
--   and runs the Newton's method on the power series.
--   The smallest ratio of the coefficients of the series will give
--   the convergence radius and the location for the nearest singularity.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree of the series,
  --   and for the precision.  Launches the tests.

    prc : character;

  begin
    new_line;
    put_line("Testing the Fabry ratio theorem ...");
    new_line;
    put_line("MENU for the working precision :");
    put_line("  1. standard double precision");
    put_line("  2. double double precision");
    put_line("  3. triple double precision");
    put_line("  4. quad double precision");
    put_line("  5. penta double precision");
    put_line("  6. octo double precision");
    put_line("  7. deca double precision");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(prc,"1234567");
    case prc is
      when '1' => Test_Standard_Fabry.Main;
      when '2' => Test_DoblDobl_Fabry.Main;
      when '3' => Test_TripDobl_Fabry.Main;
      when '4' => Test_QuadDobl_Fabry.Main;
      when '5' => Test_PentDobl_Fabry.Main;
      when '6' => Test_OctoDobl_Fabry.Main;
      when '7' => Test_DecaDobl_Fabry.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_fabry;
