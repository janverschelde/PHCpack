with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Double_Real_Powered_Series;
with Real_Powered_Series_IO;

package body Test_Real_Powered_Series is

  procedure Write ( cff : in Standard_Complex_Vectors.Vector;
                    pwt : in Standard_Floating_Vectors.Vector;
                    t : in character := 't' ) is

  -- DESCRIPTION :
  --   Plain output, for testing purposes.

  begin
    put(cff(0)); new_line;
    for k in pwt'range loop
      put(cff(k)); put(" * "); put(t); put("**");
      put(pwt(k),1,14,3); new_line;
    end loop;
  end Write;

  procedure Random_Series
              ( size : in integer32;
                cff : out Standard_Complex_Vectors.Vector;
                pwt : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Generates a random series with positive powers of the
  --   given size and sorts the powers.

  begin
    cff := Standard_Random_Vectors.Random_Vector(0,size);
    pwt := Standard_Random_Vectors.Random_Vector(1,size);
    for i in pwt'range loop
      pwt(i) := abs(pwt(i));
    end loop;
    Double_Real_Powered_Series.sort(pwt,cff);
  end Random_Series;

  procedure Test_String_Series ( size : in integer32 ) is

    cff : Standard_Complex_Vectors.Vector(0..size);
    pwt : Standard_Floating_Vectors.Vector(1..size);

  begin
    Random_Series(size,cff,pwt);
    put_line("a random series : "); Write(cff,pwt);
    declare
      s : constant string := Real_Powered_Series_IO.to_string(cff,pwt);
    begin
      new_line;
      put_line("the string format of the series :");
      put_line(s);
    end;
  end Test_String_Series;

  procedure Test_Random_Series ( size : in integer32 ) is

    cff : Standard_Complex_Vectors.Vector(0..size);
    pwt : Standard_Floating_Vectors.Vector(1..size);
    ans : character;

  begin
    Random_Series(size,cff,pwt);
    put_line("a random series : "); Write(cff,pwt);
    new_line;
    put("Write in new_line format ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y'
     then Real_Powered_Series_IO.put_line(cff,pwt);
     else Real_Powered_Series_IO.put(cff,pwt); new_line;
    end if;
  end Test_Random_Series;

  procedure Main is

    deg : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the size of the series : "); get(deg);
    new_line;
    put("-> generating a random series of size "); put(deg,1);
    put_line(" ...");
    new_line;
    put("Test output to string ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_String_Series(deg);
     else Test_Random_Series(deg);
    end if;
  end Main;

end Test_Real_Powered_Series;
