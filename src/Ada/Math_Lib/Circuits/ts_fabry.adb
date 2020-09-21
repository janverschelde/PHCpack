with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Standard_Fabry;
with Test_DoblDobl_Fabry;
with Test_QuadDobl_Fabry;

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
    put_line("Developing series starting at a regular solution ...");
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(prc,"012");
    case prc is
      when '0' => Test_Standard_Fabry.Main;
      when '1' => Test_DoblDobl_Fabry.Main;
      when '2' => Test_QuadDobl_Fabry.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_fabry;
