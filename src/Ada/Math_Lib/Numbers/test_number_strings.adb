with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Strings_and_Numbers;                use Strings_and_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;

package body Test_Number_Strings is

  procedure Standard_Test is

    use Standard_Complex_Numbers;

    f : double_float := 0.0;
    c : Complex_Number;
    s : string(1..22) := (1..22 => ' ');

  begin
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
      s2 : string(1..46);
      dp : natural32 := 0;
      s3 : string(1..46);
      c2,c3 : Complex_Number;
      last : integer;
    begin
      put("your complex number : "); put(c); new_line;
      put("     -> as a string : "); put_line(s);
      put(s2,c);
      put("->  as basic string : "); put_line(s2);
      put("Give number of decimal places : "); get(dp);
      put(s3,c,dp);
      put("-> shortened format : "); put_line(s2);
      put_line("parsing the two strings ...");
      get(s2,c2,last);
      put("-> basic format : "); put(c2); new_line;
      get(s3,c3,last);
      put("-> short format : "); put(c3); new_line;
    end;
    skip_line;
    new_line;
    put_line("Reading a string for a complex number ...");
    declare
      s : constant string := Read_String;
      cs : Complex_Number;
      last : integer;
    begin
      new_line;
      put_line("The string read : "); put(s);
      get(s,cs,last);
      new_line;
      put_line("The complex number : "); put(cs); new_line;
    end;
  end Standard_Test;

  procedure Multprec_Test is

    use Multprec_Complex_Numbers;

    f : Floating_Number;
    c : Complex_Number;

  begin
    put("Give a float f : "); get(f);
    put(" -> your float : "); put(f); new_line;
    declare
      cs : constant natural32 := Character_Size(f);
      sf : string(1..integer(cs));
      f2 : Floating_Number;
      last : integer;
    begin
      put("Number of characters in f : "); put(cs,1); new_line;
      put(sf,f);
      put("-> f in string : "); put_line(sf);
      get(sf,f2,last);
      put("-> after parse : "); put(f2); new_line;
    end;
    new_line;
    put("Give a complex number c : "); get(c);
    put(" -> your complex_number : "); put(c); new_line;
    declare
      cs : constant natural32 := Character_Size(c);
      sf : string(1..integer(cs));
      c2 : Complex_Number;
      last : integer;
    begin
      put("Number of characters in c : "); put(cs,1); new_line;
      put(sf,c);
      put("-> f in string : "); put_line(sf);
      get(sf,c2,last);
      put("-> after parse : "); put(c2); new_line;
    end;
  end Multprec_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for selecting precision in string parsing and writing :");
    put_line("  1. standard double floating-point precision; or");
    put_line("  2. arbitrary multiprecision floating-point numbers.");
    put("Type 1 or 2 to select the precision : ");
    Ask_Alternative(ans,"12");
    new_line;
    case ans is
      when '1' => Standard_Test;
      when '2' => Multprec_Test;
      when others => null;
    end case;
  end Main;

end Test_Number_Strings;
