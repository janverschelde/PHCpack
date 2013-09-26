with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Strings_and_Numbers;                use Strings_and_Numbers;

procedure ts_strnum is

-- DESCRIPTION :
--   Test on converting complex floating-point numbers to strings
--   and parsing strings back into complex floating-point numbers.

  procedure Main is

    f : double_float := 0.0;
    c : Complex_Number;
    s : string(1..21) := (1..21 => '#');

  begin
    new_line;
    put("Give a double float : "); get(f);
    put("      -> your float : "); put(f); new_line;
    put(s,f);
    put_line("the float written to a string : "); 
    put_line(s);
    declare
      s : constant string := Convert(f);
    begin
      put("your double float : "); put(f); new_line;
      put("   -> as a string : "); put_line(s);
    end;
    new_line;
    put("Give a complex number : "); get(c);
    declare
      s : constant string := Signed_Coefficient(c);
    begin
      put("your complex number : "); put(c); new_line;
      put("     -> as a string : "); put_line(s);
    end;
  end Main;

begin
  Main;
end ts_strnum;
