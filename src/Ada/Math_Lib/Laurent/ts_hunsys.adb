with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;          use Communications_with_User;
with Test_Leading_Powers;
with Test_Leading_Terms;

procedure ts_hunsys is

-- DESCRIPTION :
--   Calls the main test on leading powers of a linear system
--   of real power series.

  ans : character;

begin
  new_line;
  put("Focus test only on leading powers ? (y/n) "); Ask_Yes_or_No(ans);
  if ans ='y'
   then Test_Leading_Powers.main;
   else Test_Leading_Terms.main;
  end if;
end ts_hunsys;
