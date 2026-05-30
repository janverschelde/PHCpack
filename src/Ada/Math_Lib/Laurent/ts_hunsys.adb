with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;
with Test_Leading_Powers;
with Test_Leading_Terms;
with Test_Linear_Series_Solver;

procedure ts_hunsys is

-- DESCRIPTION :
--   Calls the main testers on solving a linear system of real power series.

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test linear system of real power series :");
    put_line("  0. test only leading powers");
    put_line("  1. test system with leading terms in coefficient matrix");
    put_line("  2. test general linear systems of real power series");
    put("Type 0, 1, or 2 to select a test : ");
    Communications_with_User.Ask_Alternative(ans,"012");
    case ans is
      when '0' => Test_Leading_Powers.main;
      when '1' => Test_Leading_Terms.main;
      when '2' => Test_Linear_Series_Solver.main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hunsys;
