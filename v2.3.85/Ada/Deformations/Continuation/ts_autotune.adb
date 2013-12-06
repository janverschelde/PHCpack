with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Continuation_Parameters;
with Continuation_Parameters_io;

procedure ts_autotune is

-- DESCRIPTION :
--   There are two master parameters to set the continuation parameters:
--   (1) the level of difficulty of the solution path;
--   (2) the number of decimal places in the working precision.
--   The higher the level, the smaller the step size.
--   The higher the number of decimal places, the smaller the tolerances.

  procedure Main is

    level : natural32 := 0;
    nbdgts : natural32 := 16;
    ans : character;

  begin
    new_line;
    put_line("Tuning the continuation parameters ...");
    loop
      new_line;
      put_line("Current values of the two master parameters :");
      put("  difficulty level : "); put(level,1); new_line;
      put("  number of digits : "); put(nbdgts,1); new_line;
      new_line;
      Continuation_Parameters.Tune(level,nbdgts);
      Continuation_Parameters_io.put;
      new_line;
      put("Change the current settings ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      new_line;
      put("Give new difficulty level : "); get(level);
      put("Give new number of digits : "); get(nbdgts);
    end loop;
  end Main;

begin
  Main;
end ts_autotune;
