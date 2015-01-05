with unchecked_deallocation;
with integer_io;                         use integer_io;
with File_Scanning;                      use File_Scanning;

package body String_Splitters is

  function Read_String return string is

    temp : string(1..256);
    cnt : natural;

  begin
    put("Give a string of characters : ");
    get_line(temp,cnt);
    return temp(1..cnt);
  end Read_String;

  procedure Append ( s : in out Link_to_String; t : in string ) is
  begin
    if s = null then
      s := new string'(t);
    else
      declare
        ns : constant string := s.all & t;
      begin
        Clear(s);
        s := new string'(ns);
      end;
    end if;
  end Append;

  function Read_till_Delimiter
             ( file : in file_type; d : in character )
             return Link_to_String is

    res : Link_to_String := null;
    buffer : string(1..256) := (1..256 => ' ');
    c : character;
    ind : natural := 0;

  begin
    while not End_of_File(file) loop
      get(file,c);
      exit when (c = d);
      ind := ind + 1;
      if ind > buffer'last then
        Append(res,buffer);
        ind := 1;
      end if;
      buffer(ind) := c; 
    end loop;
    if c = d then
      ind := ind + 1;
      buffer(ind) := c;
    end if;
    Append(res,buffer(1..ind));
    return res;
  end Read_till_Delimiter;

  function Read_till_Delimiter
             ( file : in file_type; n : in natural; d : in character )
             return Array_of_Strings is

    res : Array_of_Strings(1..n);

  begin
    for i in 1..n loop
      if End_of_File(file)
       then res(i) := null;
       else res(i) := Read_till_Delimiter(file,d);
      end if;
    end loop;
    return res;
  end Read_till_Delimiter;

  procedure get ( n,m : out natural; p : out Link_to_Array_of_Strings ) is
  begin
    get(standard_input,n,m,p);
  end get;

  procedure get ( file : in file_type;
                  n,m : out natural; p : out Link_to_Array_of_Strings ) is
  begin
    get(file,n);
    m := Scan_Line_for_Number(file);
    if m = 0
     then m := n;
    end if;
    declare
      s : Array_of_Strings(1..n);
    begin
      s := Read_till_Delimiter(file,n,';');
      p := new Array_of_Strings'(s);
    end;
  end get;

  function Count_Delimiters ( s : string; d : character ) return natural is

    res : natural := 0;

  begin
    for i in s'range loop
      if s(i) = d
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Count_Delimiters;

  function Split ( n : natural; s : string; d : character )
                 return Array_of_Strings is

    res : Array_of_Strings(1..n);
    ind : integer := 0;
    buf : string(s'range);
    cnt : integer := buf'first-1;

  begin
    for i in s'range loop
      cnt := cnt + 1;
      buf(cnt) := s(i);
      if s(i) = d then
        ind := ind + 1;
        res(ind) := new string'(buf(buf'first..cnt));
        cnt := buf'first-1;
      end if;
    end loop;
    return res;
  end Split;

-- DESTRUCTORS :

  procedure Clear ( s : in out Link_to_String ) is

    procedure free is new unchecked_deallocation(string,Link_to_String);

  begin
    if s /= null
     then free(s);
    end if;
  end Clear;

  procedure Clear ( s : in out Array_of_Strings ) is
  begin
    for i in s'range loop
      Clear(s(i));
    end loop;
  end Clear;

  procedure Clear ( s : in out Link_to_Array_of_Strings ) is

    procedure free is
      new unchecked_deallocation(Array_of_Strings,Link_to_Array_of_Strings);

  begin
    if s /= null then
      for i in s'range loop
        Clear(s(i));
      end loop;
      free(s);
    end if;
  end Clear;

end String_Splitters;
