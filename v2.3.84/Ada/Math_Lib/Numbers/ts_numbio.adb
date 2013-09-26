with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Parse_Numbers;              use Standard_Parse_Numbers;

procedure ts_numbio is

-- DESCRIPTION :
--   Tests on the interactive reading of numbers with
--   appropriate exception handling.

  function Read_String ( file : in file_type ) return string is

  -- DESCRIPTION :
  --   Reads a string from file, suppressing exceptions.

    tmp : string(1..256);
    cnt : natural;

  begin
    get_line(file,tmp,cnt);
    return tmp(1..cnt);
  exception
    when others
      => put_line("reading of characters failed"); return "";
  end Read_String;

  procedure Write_Error_Message ( s : in string ) is

  -- DESCRIPTION :
  --   Writes the error message showing the string.

  begin
    new_line; put("'"); put(s); put("'");
    put_line(" contains no natural number");
  end Write_Error_Message;

  procedure Check_Trailing_Characters
              ( s : in string; p : in integer; fail : out boolean ) is

  -- DESCRIPTION :
  --   Returns fail as false if s(p..s'last) contains only spaces,
  --   returns fails as true otherwise.

  begin
    for i in p..s'last loop
      fail := not (s(i) = ' '); 
      exit when fail;
    end loop;
    if fail then
      new_line;
      put("The trailing characters '"); put(s(p..s'last));
      put("' in '"); put(s); put_line("' are not spaces");
    end if;
  end Check_Trailing_Characters;

  procedure Read_Natural ( n : out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Reads a string from file, parsed the string into a natural number
  --   and writes error messages if the parsing fails.
  --   If not fail, then n contains the natural number given by the user.

  begin
    fail := true; n := 0;
    put("Give a natural number : ");
    declare
      s : constant string := Read_String(standard_input);
      p : integer := s'first;
      sign : character;
      i : integer32;
      ni : natural32;
    begin
      Parse(s,p,i,ni,sign);
     -- put("parsed "); put(ni,1); put_line(" characters");
      if ni = 0 then
        Write_Error_Message(s); fail := true;
      elsif i < 0 then
        put("the parsed number "); put(i,1); put_line(" is negative");
        fail := true;
      else
        fail := false;
        if p <= s'last
         then Check_Trailing_Characters(s,p,fail);
        end if;
        if not fail
         then n := natural32(i);
        end if;
      end if;
    exception
      when others => Write_Error_Message(s);
    end;
  end Read_Natural;

  procedure Main is

    n : natural32;
    fail : boolean;

  begin
    new_line;
    put_line("Reading a natural number ...");
    Read_Natural(n,fail);
    if fail
     then put_line("The reading of a natural number failed.");
     else put("Read the natural number "); put(n,1); put_line(".");
    end if;
  end Main;

begin
  Main;
end ts_numbio;
