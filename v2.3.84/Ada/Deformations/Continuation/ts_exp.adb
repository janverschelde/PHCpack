with text_io;                             use text_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;     use Standard_Mathematical_Functions;

procedure ts_exp is

  a,b : double_float := 0.0;
  ans : character;

begin
  new_line;
  put_line("Calculating a**b, for a,b floating points");
  new_line;
  loop
    put("Give a : "); get(a);
    loop
      put("Give b : "); get(b);
      put(a); put("**"); put(b);
      put(" : "); put(a**b); new_line;
      put("Do you want other exponents of a ? (y/n) ");
      get(ans);
      exit when (ans /= 'y');
    end loop;
    put("Do you want other values for a ? (y/n) ");
    get(ans);
    exit when (ans /= 'y');
  end loop;
end ts_exp;
