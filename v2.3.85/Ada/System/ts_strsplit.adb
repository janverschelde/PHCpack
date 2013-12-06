with text_io,integer_io;                 use text_io,integer_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;                   use String_Splitters;

procedure ts_strsplit is

-- DESCRIPTION :
--   Test the operations on strings.

  procedure Test_Read_one_String is
  begin
    new_line;
    put_line("Reading one string ...");
    declare
      s : constant string := Read_String;
    begin
      put("read '"); put(s); put_line("'");
    end;
  end Test_Read_one_String;

  procedure Test_Read_till_Delimiter is
  begin
    new_line;
    put_line("Give a string, terminated by semicolon ...");
    declare
      p : constant Link_to_String := Read_till_Delimiter(standard_input,';');
    begin
      put("your string is '"); put(p.all); put_line("'");
    end;
  end Test_Read_till_Delimiter;

  procedure Test_Read_Sequence is

  -- DESCRIPTION :
  --   Prompts the user for the name of a file, reads n from that file,
  --   and then a sequence of n polynomials terminated by a semicolon.

    file : file_type;
    n : natural;

  begin
    new_line;
    put_line("Reading a file name for a polynomial system...");
    Read_Name_and_Open_File(file);
    get(file,n);
    new_line;
    put("Reading "); put(n,1); put_line(" strings from file...");
    declare
      p : constant Array_of_Strings(1..n) := Read_till_Delimiter(file,n,';');
    begin
      for i in 1..n loop
        put("p["); put(i,1); put("] : "); put(p(i).all); new_line;
      end loop;
    end;
  end Test_Read_Sequence;

  procedure Test_Read_System is

  -- DESCRIPTION :
  --   The system read from file can have a different number of
  --   equations than a number of variables.

    file : file_type;
    p : Link_to_Array_of_Strings;
    n,m : natural;

  begin
    new_line;
    put_line("Reading a file name for a polynomial system...");
    Read_Name_and_Open_File(file);
    get(file,n,m,p);
    new_line;
    put("number of equations : "); put(n,1); new_line;
    put("number of variables : "); put(m,1); new_line;
    for i in 1..n loop
      put("p["); put(i,1); put("] : "); put(p(i).all); new_line;
    end loop;
  end Test_Read_System;

  procedure Main is
 
    ans : character;
   
  begin
    new_line;
    put_line("MENU for testing string reading : ");
    put_line("  1. reading one string from standard input on a line;");
    put_line("  2. reading one string from standard input till delimiter;");
    put_line("  3. reading a sequence of polynomials from file.");
    put_line("  4. reading a system of polynomials from file.");
    put("Type 1, 2, 3, or 4 to choose : "); Ask_Alternative(ans,"1234");
    case ans is 
      when '1' => Test_Read_one_String;
      when '2' => Test_Read_till_Delimiter;
      when '3' => Test_Read_Sequence;
      when others => Test_Read_System;
    end case;
  end Main;

begin
  Main;
end ts_strsplit;
