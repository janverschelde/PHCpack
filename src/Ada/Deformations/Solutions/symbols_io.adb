with Symbol_Table_io;

package body Symbols_io is

  procedure Skip_Symbol ( file : in file_type ) is

    c : character;

  begin
    loop
      get(file,c);
      exit when (c = ':');
    end loop;
  end Skip_Symbol;

  function Read_Symbol ( file : file_type ) return Symbol is

    sb : Symbol;
    c : character;

  begin
    for i in sb'range loop
      sb(i) := ' ';
    end loop;
    loop       -- skip the spaces
      get(file,c);
      exit when ((c /= ' ') and (c /= ASCII.CR));
    end loop;
    sb(1) := c;
    for i in sb'first+1..sb'last loop
      get(file,c);
      exit when c = ' ';
      sb(i) := c;
    end loop;
    return sb;
  end Read_Symbol;    

  function Get_Symbol ( file : file_type ) return natural32 is

    sb : constant Symbol := Read_Symbol(file);

  begin
    return Symbol_Table.get(sb);
  end Get_Symbol;

  function Get_Symbol ( file : file_type;
                        s : Array_of_Symbols ) return natural32 is

    sb : constant Symbol := Read_Symbol(file);

  begin
    for i in s'range loop
      if Equal(s(i),sb)
       then return natural32(i);
      end if;
    end loop;
    return 0;
  end Get_Symbol;

  procedure put_symbol ( file : in file_type; i : in natural32 ) is

    sb : constant Symbol := Get(i);

  begin
    Symbol_Table_io.put(file,sb);
  end put_symbol;

end Symbols_io;
