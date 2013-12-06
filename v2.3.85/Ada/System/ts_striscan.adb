with text_io,integer_io;                 use text_io,integer_io;
with String_Splitters;                   use String_Splitters;
with String_Parsing;

procedure ts_striscan is

-- DESCRIPTION :
--   Interactive test on scanning for banners in a string.

  procedure Main is

    s : constant string := Read_String;
    b : constant string := Read_String;
    i : constant integer := String_Parsing.Scan(s,b);

  begin
    if i < 0 then
      put_line(b & " does not occur in " & s);
    else
      put_line(b & " occurs in " & s);
      put("end position : "); put(i,1); new_line;
      put("sub string : "); put_line(s(i-b'last+1..i));
    end if;
  end Main;

begin
  new_line;
  put_line("Scanning string for a banner.");
  new_line;
  put_line("Reading a string and a banner...");
  Main;
end ts_striscan;
