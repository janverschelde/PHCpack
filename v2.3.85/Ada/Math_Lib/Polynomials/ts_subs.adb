with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Substitutors;     use Standard_Complex_Substitutors;

procedure ts_subs is

-- DESCRIPTION :
--   Substitutes one variable by a random number.

  procedure Interactive_Substitute
              ( file : in file_type; p : in Poly_Sys ) is

    ans : character;
    r : Complex_Number;

  begin
    for i in 1..Symbol_Table.Number loop
      put("Variable "); put(i,1); put(" is ");
      Symbol_Table_io.put(Symbol_Table.get(i)); new_line;
      put("Substitute variable "); put(i,1); put(" ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Symbol_Table.Remove(i);
        r := Random1;
        declare
          rp : constant Poly_Sys(p'range) := Substitute(integer32(i),r,p);
        begin
          put_line(file,rp);
          new_line;
          put_line("See the output file for results...");
          new_line;
        end;
      end if;
      exit when ans = 'y';
    end loop;
  end Interactive_Substitute;

  procedure Main is

    lp : Link_to_Poly_Sys;
    file : file_type;

  begin
    new_line;
    put_line("Substitution of one variable by a random constant.");
    new_line;
    put_line("Reading the input system...");
    get(lp);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    Interactive_Substitute(file,lp.all);
  end Main;

begin
  Main;
end ts_subs;
