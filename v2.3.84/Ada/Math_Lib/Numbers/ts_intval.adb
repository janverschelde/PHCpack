with text_io;                          use text_io;
with Communications_with_User;         use Communications_with_User;
with Standard_Interval_Numbers;        use Standard_Interval_Numbers;
with Standard_Interval_Numbers_io;     use Standard_Interval_Numbers_io;
with Standard_Interval_Ring;
with Standard_Interval_Ring_io;
with Standard_Interval_Ring.FField;

procedure ts_intval is

-- DESCRIPTION :
--   Interactive development of interval arithmetic.

  procedure Test_IO is

  -- DESCRIPTION :
  --   Tests basic input/output and thus also create.

    x : Interval;

  begin
    put("Give interval : "); get(x);
    put("Your interval : "); put(x); new_line;
  end Test_IO;

  procedure Test_Arithmetic is

  -- DESCRIPTION :
  --   Prompts the user for two intervals and shows the
  --   result of all four basic operations.

    x,y : Interval;

  begin
    put("Give interval x : "); get(x); skip_line;
    put("x = "); put(x); new_line;
    put("Give interval y : "); get(y); skip_line;
    put("y = "); put(y); new_line;
    put("x + y = "); put(x+y); new_line;
    put("x - y = "); put(x-y); new_line;
    put("x * y = "); put(x*y); new_line;
    put("x / y = "); put(x/y); new_line;
  end Test_Arithmetic;

  procedure Main is

    ans : character;
 
  begin
    new_line;
    put_line("Test on interval arithmetic ...");
    new_line;
    loop
      put_line("Choose one of the following operations : ");
      put_line("  0. exit this program;");
      put_line("  1. test basic input/output;");
      put_line("  2. test interval arithmetic operations.");
      put("Type 0, 1, or 2 to select : "); Ask_Alternative(ans,"012");
      exit when ans = '0';
      new_line;
      if ans = '1' then
        Test_IO;
      else
        Test_Arithmetic;
      end if;
    end loop;
  end Main;

begin
  Main;
end ts_intval;
