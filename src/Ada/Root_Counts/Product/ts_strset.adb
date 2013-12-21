with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Random_Product_Start_Systems;
with Set_Structure;
with Set_Structure_io;
with Set_Structure_Strings;
with Supporting_Set_Structure;

procedure ts_strset is

-- DESCRIPTION :
--   Interactive test on string representations of set structures.

  procedure Strings_and_Parsing is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, makes a set structure,
  --   writes its string representation and then parses this string. 

    p : Link_to_Poly_Sys;

  begin
    new_line;
    get(p);
    new_line;
    Random_Product_Start_Systems.Build_Set_Structure(p.all);
    put_line("a supporting set structure :");
    Set_Structure_io.put;
    declare
      s : constant string := Set_Structure_Strings.to_string;
    begin
      put_line("its string representation :");
      put_line(s);
      put_line("clearing the set structure ...");
      Set_Structure.Clear;
      put_line("parsing the set structure string ...");
      Set_Structure_Strings.Parse(s);
      put_line("after parsing the set structure string :");
      Set_Structure_io.put;
    end;
  end Strings_and_Parsing;

  procedure Check_if_Supporting is

    p : Link_to_Poly_Sys;

  begin
    new_line;
    get(p);
    new_line;
    put_line("Reading a set structure string ...");
    declare
      s : constant string := String_Splitters.Read_String;
      ans : character;
      otp : boolean;
    begin
      put_line("Your set structure : ");
      put_line(s);
      Set_Structure_Strings.Parse(s);
      put_line("after parsing the set structure string :");
      Set_Structure_io.put;
      new_line;
      put("Do you want intermediate output ? (y/n) ");
      Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      new_line;
      if Supporting_Set_Structure.Is_Supporting(p.all,otp)
       then put_line("The set structure supports the system.");
       else put_line("The set structure does not support the system.");
      end if;
    end;
  end Check_if_Supporting;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU on strings for supporting set structures :");
    put_line("  1. test on writing to strings and parsing from strings;");
    put_line("  2. check whether set structure is supporting.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Strings_and_Parsing;
     else Check_if_Supporting;
    end if;
  end Main;

begin
  Main;
end ts_strset;
