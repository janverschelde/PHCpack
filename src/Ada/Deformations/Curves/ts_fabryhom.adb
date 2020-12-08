with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Newton_Fabry_on_Homotopy;

procedure ts_fabryhom is

-- DESCRIPTION :
--   Tests the Newton-Fabry convergence radius computation
--   for artificial or natural-parameter homotopies.

  ans : character;

begin
  new_line;
  put("Generate a test homotopy ? (y/n) "); Ask_Yes_or_No(ans);
  if ans /= 'y'
   then Newton_Fabry_on_Homotopy.Main;
   else Newton_Fabry_on_Homotopy.Test;
  end if;
end ts_fabryhom;
