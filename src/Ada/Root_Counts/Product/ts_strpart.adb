with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Symbol_Table,Symbol_Table_io;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Sets_of_Unknowns_io;                use Sets_of_Unknowns_io;
with Sets_of_Unknowns_Strings;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;
 use Partitions_of_Sets_of_Unknowns_io;
with Partitions_of_Sets_Strings;

procedure ts_strpart is

-- DESCRIPTION :
--   Test on strings and partitions of sets of unknowns.

  procedure Initialize_Symbol_Table ( n : in natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for n symbols and initializes the symbol table.

    use Symbol_Table;
    sbs : Array_of_Symbols(1..integer32(n));

  begin
    for i in sbs'range loop
      put("Give symbol "); put(i,1); put(" : ");
      Symbol_Table_io.get(sbs(i)); skip_line;
    end loop;
    put("The symbols : ");
    for i in sbs'range loop
      put(" "); Symbol_Table_io.put(sbs(i));
    end loop;
    new_line;
    Symbol_Table.Init(sbs);
  end Initialize_Symbol_Table;

  procedure Test_Strings_of_Sets is

  -- DESCRIPTION :
  --   Prompts the user for a number of variables and then
  --   builds the set adding one variable after the other,
  --   each time printing the current set.

    n : natural32 := 0;
    s : Set;

  begin
    put("Give the number of variables : "); get(n);
    s := Create(n);
    for i in 1..n loop
      put("Set before adding variable "); put(i,1); put(" : "); put(s);
      new_line;
      put("-> the string representation : "); 
      put_line(Sets_of_Unknowns_Strings.to_string(s));
      Add(s,i);
    end loop;
    put("The final set : "); put(s); new_line;
    put("-> the string representation : "); 
    put_line(Sets_of_Unknowns_Strings.to_string(s));
  end Test_Strings_of_Sets;

  procedure Test_Parse_Strings_into_Sets is

  -- DESCRIPTION :
  --   Prompts the user for a string and then parses the string
  --   into a set of unknowns.

    n : natural32 := Symbol_Table.Number;
    ans : character;

  begin
    new_line;
    if n = 0 then
      put("Give the number of variables : "); get(n); skip_line;
      new_line;
      put("Initialize the symbol table ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Initialize_Symbol_Table(n);
      end if;
      new_line;
    end if;
    put_line("Reading a string with a set of unknowns ...");
    declare
      s : constant string := String_Splitters.Read_String;
      r : Set;
    begin
      put_line("The string read : " & s);
      r := Sets_of_Unknowns_Strings.Parse(s,n);
      put("The set parsed : "); put(r); new_line;
    end;
  end Test_Parse_Strings_into_Sets;

  procedure Show_all_Partitions ( s : in Set ) is

  -- DESCRIPTION :
  --   Writes all partitions of a set of unknowns to screen.

    cnt : natural32 := 0;

    procedure Write ( p : in Partition; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      put("partition "); put(cnt,1); put(" : "); put(p); new_line;
      put("-> its string representation : ");
      put_line(Partitions_of_Sets_Strings.to_string(p));
      continue := true;
    end Write;
    procedure Write_Partitions is new Generate_Partitions(Write);

  begin
    Write_Partitions(s);
  end Show_all_Partitions;

  procedure Test_Strings_of_Partitions is

  -- DESCRIPTION :
  --   Prompts the user for a number of variables in the set
  --   and then writes all partitions of that set.

    n : natural32 := 0;
    s : Set;

  begin
    put("Give the number of variables : "); get(n);
    s := Create(n);
    for i in 1..n loop
      Add(s,i);
    end loop;
    put("The set : "); put(s); new_line;
    new_line;
    Show_all_Partitions(s);
  end Test_Strings_of_Partitions;

  procedure Test_Parse_Strings_into_Partitions is

    n : natural32 := Symbol_Table.Number;
    ans : character;

  begin
    new_line;
    if n = 0 then
      put("Give the number of variables : "); get(n); skip_line;
      new_line;
      put("Initialize the symbol table ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Initialize_Symbol_Table(n);
      end if;
      new_line;
    end if;
    put_line("Reading a string with a partition of a set of unknowns ...");
    declare
      s : constant string := String_Splitters.Read_String;
      r : constant Partition := Partitions_of_Sets_Strings.Parse(s,n);
    begin
      put_line("The string read : " & s);
      put("The partition parsed : "); put(r); new_line;
    end;
  end Test_Parse_Strings_into_Partitions;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test string operations on sets and partitions :");
    put_line("  1. write string representations for sets;");
    put_line("  2. parse string into a set of unknowns;");
    put_line("  3. write string representations for partitions;");
    put_line("  4. parse string into a partition of a set of unknowns.");
    put("Type 1, 2, 3, or 4 to make your choice : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Test_Strings_of_Sets;
      when '2' => Test_Parse_Strings_into_Sets;
      when '3' => Test_Strings_of_Partitions;
      when '4' => Test_Parse_Strings_into_Partitions;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_strpart;
