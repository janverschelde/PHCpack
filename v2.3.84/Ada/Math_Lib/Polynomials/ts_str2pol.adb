with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Characters_and_Numbers;
with Parse_Strings_to_Polynomials;       use Parse_Strings_to_Polynomials;
with Symbol_Table;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Strings;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions; 

procedure ts_str2pol is

-- DESCRIPTION :
--   Development of parsing a sequence of strings into polynomials,
--   with error messages written to strings as well.

  function new_line_string return string is

  -- DESCRIPTION :
  --   Returns the new line symbol as a string.

    res : string(1..1);

  begin
    res(1) := ASCII.LF;
    return res;
  end new_line_string;

  function to_string ( sb : Symbol_Table.Symbol ) return string is

  -- DESCRIPTION :
  --   Returns the symbol (without trailing spaces) and preceded
  --   with one space.

    res : string(1..80);
    cnt : natural := 1;

  begin
    res(1) := ' ';
    for i in sb'range loop
      exit when (sb(i) = ' ');
      cnt := cnt + 1;
      res(cnt) := sb(i);
    end loop;
    return res(1..cnt);
  end to_string;

  procedure Append_Symbols
              ( s : out String_Splitters.Link_to_String ) is

  -- DESCRIPTION :
  --   Appends the content of the symbol table to the string s.

  begin
    for i in 1..Symbol_Table.Number loop
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.get(i);
        s_sb : constant string := to_string(sb);
      begin
        String_Splitters.Append(s,s_sb);
      end;
    end loop;
  end Append_Symbols;

  procedure Parse_Polynomials
              ( n : in natural32;
                a : in String_Splitters.Link_to_Array_of_Strings;
                p : in Link_to_Laur_Sys; fail : out integer32;
                emsg : out String_Splitters.Link_to_String ) is

  -- DESCRIPTION :
  --   Parses the strings given in a into polynomials.

  -- REQUIRED : a'range = p'range.

  -- ON ENTRY :
  --   n        number of variables;
  --   a        input strings are supposed to contain polynomials.

  -- ON RETURN :
  --   p        parsed polynomials if fail differs from zero;
  --   fail     index of the first string in a where the parsing
  --            raised an exception;
  --   emsg     error message if fail index > 0.

    procedure Write_Error ( index : in integer32 ) is

    -- DESCRIPTION :
    --   Writes a generic error message when the error occurs
    --   in the string with the given index.

      s_index : constant string
              := Characters_and_Numbers.nConvert(natural32(index));

    begin
      String_Splitters.Append(emsg,"Error raised in string ");
      String_Splitters.Append(emsg,s_index);
      String_Splitters.Append(emsg," the wrong input is on next line");
      String_Splitters.Append(emsg,new_line_string);
      String_Splitters.Append(emsg,a(integer(index)).all);
      String_Splitters.Append(emsg,new_line_string);
    end Write_Error;

  begin
    fail := 0; 
    emsg := new string'("");
    for i in a'range loop
      exit when (fail > 0);
      declare
        q : Poly;
        s : constant string := a(i).all;
        s_index : constant string 
                := Characters_and_Numbers.nConvert(natural32(i));
      begin
        q := Standard_Complex_Laur_Strings.Parse(n,s);
        p(integer32(i)) := q;
      exception
        when ILLEGAL_CHARACTER =>
          fail := integer32(i); Write_Error(integer32(i));
          String_Splitters.Append(emsg,
            "An invalid character was found at polynomial ");
          String_Splitters.Append(emsg,s_index);
          String_Splitters.Append(emsg,".");
          String_Splitters.Append(emsg,new_line_string);
        when ILLEGAL_OPERATION =>
          fail := integer32(i); Write_Error(integer32(i));
          String_Splitters.Append(emsg,
            "An invalid operation was found at polynomial ");
          String_Splitters.Append(emsg,s_index);
          String_Splitters.Append(emsg,".");
          String_Splitters.Append(emsg,new_line_string);
        when OVERFLOW_OF_UNKNOWNS =>
          fail := integer32(i); Write_Error(integer32(i));
          String_Splitters.Append(emsg,
            "Too many unknowns at polynomial ");
          String_Splitters.Append(emsg,s_index);
          String_Splitters.Append(emsg,".");
          String_Splitters.Append(emsg,new_line_string);
          String_Splitters.Append(emsg,"The current symbol table : ");
          Append_Symbols(emsg);
          String_Splitters.Append(emsg,new_line_string);
        when BAD_BRACKET =>
          fail := integer32(i); Write_Error(integer32(i));
          String_Splitters.Append(emsg,"Misplaced bracket at polynomial ");
          String_Splitters.Append(emsg,s_index);
          String_Splitters.Append(emsg,".");
          String_Splitters.Append(emsg,new_line_string);
        when others =>
          fail := integer32(i); Write_Error(integer32(i));
      end;
    end loop;
  end Parse_Polynomials;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a file name
  --   and then parses the strings from that file.

    infile,outfile : file_type;
    nq,nv : natural32 := 0;
    fail : boolean;
    las : String_Splitters.Link_to_Array_of_Strings;
    lq : Link_to_Laur_Sys;
    fail_index : integer32;
    error : String_Splitters.Link_to_String;

  begin
    new_line;
    put_line("Reading the name of the input file...");
    Read_Name_and_Open_File(infile);
    new_line;
    Read_First_Line(infile,outfile,false,nq,nv,fail);
    if fail then
      put_line("Parsing first line of input for dimensions failed.");
    else
      new_line;
      put("the number of polynomials : "); put(nq,1); new_line;
      put("the number of variables : "); put(nv,1); new_line;
      Read_Polynomials(infile,outfile,false,nq,las,fail);
      new_line;
      if fail then
        put("Failed to read "); put(nq,1);
        put_line(" strings separated by semicolons.");
      else
        put("Read "); put(nq,1); put_line(" strings from file :");
        for i in las'range loop
          put_line(las(i).all);
        end loop;
        Symbol_Table.Init(nv);
        lq := new Laur_Sys(1..integer32(nq));
        Parse_Polynomials(nv,las,lq,fail_index,error);
        new_line;
        if fail_index = 0 then
          put_line("Parsing string was okay.");
        else
          put("Error parsing string "); put(fail_index,1); put_line(".");
          new_line;
          put_line("The error message :"); put_line(error.all);
        end if;
      end if;
    end if;
  end Main;

begin
  Main;
end ts_str2pol;
