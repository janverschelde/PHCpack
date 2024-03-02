with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Greeting_Banners;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_Strings;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Black_Box_Solvers;

procedure blackbox is

-- DESCRIPTION :
--   Illustrates a basic application of the blackbox solver,
--   on a Laurent polynomial system.

  procedure Main is

  -- DESCRIPTION :
  --   Solves a Laurent polynomial system given as a string.

    s : constant string := "x*y + 2*x + 3; x^2 + y^(-2) + 4;";
    p : constant Standard_Complex_Laur_Systems.Laur_Sys
      := Standard_Complex_Laur_Strings.Parse(2,2,s);
    silent : constant boolean := false;
    rc : natural32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    put_line(Greeting_Banners.welcome);
    put_line("An example polynomial system :");
    Standard_Complex_Laur_Systems_io.put(p);
    Black_Box_Solvers.Solve(p, silent, rc, sols);
    put("The root count : "); put(rc,1); new_line;
    put_line("The solutions :");
    Standard_Complex_Solutions_io.write(standard_output,sols);
  end Main;

begin
  Main;
end blackbox;
