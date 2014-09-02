with Unchecked_Conversion;

package body Characters_and_Numbers is

  subtype short_int is integer range 0..255;

  function Integer_to_Character ( i : integer32 ) return character is

    function short_int_to_character is new
      Unchecked_Conversion(short_int,character);

  begin
    return short_int_to_character(short_int(i));
  end Integer_to_Character;

  function Character_to_Integer ( c : character ) return integer32 is

    function character_to_short_int is new
      Unchecked_Conversion(character,short_int);

  begin
    return integer32(character_to_short_int(c));
  end Character_to_Integer;

  function Convert ( c : character ) return natural32 is
  begin
    case c is
      when '0' => return 0;
      when '1' => return 1;
      when '2' => return 2;
      when '3' => return 3;
      when '4' => return 4;
      when '5' => return 5;
      when '6' => return 6;
      when '7' => return 7;
      when '8' => return 8;
      when '9' => return 9;
      when others => return 10;
    end case;
  end Convert;

  function Convert_Decimal ( n : natural32 ) return character is
  begin
    case n is
      when 0 => return '0';
      when 1 => return '1';
      when 2 => return '2';
      when 3 => return '3';
      when 4 => return '4';
      when 5 => return '5';
      when 6 => return '6';
      when 7 => return '7';
      when 8 => return '8';
      when 9 => return '9';
      when others => return '0';
    end case;
  end Convert_Decimal;

  function Convert_Hexadecimal ( c : character ) return natural32 is
  begin
     case c is
       when '0' => return 0;
       when '1' => return 1;
       when '2' => return 2;
       when '3' => return 3;
       when '4' => return 4;
       when '5' => return 5;
       when '6' => return 6;
       when '7' => return 7;
       when '8' => return 8;
       when '9' => return 9;
       when 'A' => return 10;
       when 'B' => return 11;
       when 'C' => return 12;
       when 'D' => return 13;
       when 'E' => return 14;
       when 'F' => return 15;
       when others => return 16;
     end case;
  end Convert_Hexadecimal;

  function Convert_Hexadecimal ( n : natural32 ) return character is
  begin
    if n < 10 then
      return Convert_Decimal(n);
    else
      case n is
        when 10 => return 'A';
        when 11 => return 'B';
        when 12 => return 'C';
        when 13 => return 'D';
        when 14 => return 'E';
        when 15 => return 'F';
        when others => return '0';
      end case;
    end if;
  end Convert_Hexadecimal;

  function Convert ( s : string ) return natural32 is

    acc : natural32 := 0;
    cvn : natural32;

  begin
    for i in s'range loop
      cvn := Convert(s(i));
      if cvn < 10
       then acc := acc*10 + cvn;
      end if;
    end loop;
    return acc;
  end Convert;

  function Convert ( s : string ) return natural64 is

    acc : natural64 := 0;
    cvn : natural64;

  begin
    for i in s'range loop
      cvn := natural64(Convert(s(i)));
      if cvn < 10
       then acc := acc*10 + cvn;
      end if;
    end loop;
    return acc;
  end Convert;

  function Convert ( s : string ) return integer32 is

    res : integer32 := 0;
    resnat : natural32 := 0;
    pos : integer := s'first;

  begin
    while s(pos) = ' ' loop
      pos := pos + 1;
      exit when (pos > s'last);
    end loop;
    if pos <= s'last then
      if s(pos) = '-' then
        resnat := convert(s(pos+1..s'last));
        res := -integer32(resnat);
      else
        resnat := convert(s(pos..s'last));
        res := integer32(resnat);
      end if;
    end if;
    return res;
  end Convert;

  function Convert ( s : string ) return integer64 is

    res : integer64 := 0;
    resnat : natural64 := 0;
    pos : integer := s'first;

  begin
    while s(pos) = ' ' loop
      pos := pos + 1;
      exit when (pos > s'last);
    end loop;
    if pos <= s'last then
      if s(pos) = '-' then
        resnat := convert(s(pos+1..s'last));
        res := -integer64(resnat);
      else
        resnat := convert(s(pos..s'last));
        res := integer64(resnat);
      end if;
    end if;
    return res;
  end Convert;

  function Convert ( s : string ) return double_float is

    acc : double_float := 0.0;
    cvn : natural32;

  begin
    for i in s'range loop
      cvn := Convert(s(i));
      if cvn < 10
       then acc := acc*10.0 + double_float(cvn);
      end if;
    end loop;
    return acc;
  end Convert;

  function nConvert ( n : natural32 ) return string is
  begin
    if n < 10 then
      declare
        res : String(1..1);
      begin
        res(1) := Convert_Decimal(n);
        return res;
      end;
    else
      declare
        rest : constant natural32 := n mod 10;
        head : constant natural32 := n/10;
        headstr : constant string := Convert(integer32(head));
        res : String(1..headstr'last+1);
      begin
        res(headstr'range) := headstr;
        res(res'last) := Convert_Decimal(rest);
        return res;
      end;
    end if;
  end nConvert;

  function nConvert ( n : natural64 ) return string is
  begin
    if n < 10 then
      declare
        res : String(1..1);
      begin
        res(1) := Convert_Decimal(natural32(n));
        return res;
      end;
    else
      declare
        rest : constant natural64 := n mod 10;
        head : constant natural64 := n/10;
        headstr : constant string := nConvert(head);
        res : String(1..headstr'last+1);
      begin
        res(headstr'range) := headstr;
        res(res'last) := Convert_Decimal(natural32(rest));
        return res;
      end;
    end if;
  end nConvert;

  function Convert ( i : integer32 ) return string is
  begin
    if i < 0
     then return "-" & nConvert(natural32(-i)); 
     else return nConvert(natural32(i));
    end if;
  end Convert;

  function Convert ( i : integer64 ) return string is
  begin
    if i < 0
     then return "-" & nConvert(natural64(-i)); 
     else return nConvert(natural64(i));
    end if;
  end Convert;

  procedure Skip_Spaces ( file : in file_type; c : in out character ) is
  begin
    get(file,c);
    while c = ' ' and not End_of_Line(file) loop
      get(file,c);
    end loop;
  end Skip_Spaces;

  procedure Skip_Underscores ( file : in file_type; c : in out character ) is
  begin
    get(file,c);
    while c = '_' and not End_of_Line(file) loop
      get(file,c);
    end loop;
  end Skip_Underscores;

end Characters_and_Numbers;
