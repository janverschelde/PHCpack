package body Lists_of_Symbols is

--  procedure Write_Symbols ( file : in file_type; n : in natural32 ) is
--
--    sb : Symbol;
--
--  begin
--    put(file,"[");
--    for i in 1..n loop
--      sb := Symbol_Table.Get(i);
--      for j in sb'range loop
--        exit when sb(j) = ' ';
--        put(file,sb(j));
--      end loop;
--      if i = n
--       then put_line(file," ] ,");
--       else put(file," , ");
--      end if;
--    end loop;
--  end Write_Symbols;

  procedure Create_Symbol_Table ( L : in Symbol_List ) is

    tmp : Symbol_List := L;

  begin
    Symbol_Table.Init(Length_Of(L));
    while not Is_Null(tmp) loop
      Symbol_Table.Add(Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end Create_Symbol_Table;

  function Equal ( sb : Symbol; s : string ) return boolean is
  begin
    for i in sb'range loop
      if i > s'last then
        if sb(i) = ' '
         then return true;
         else return false;
        end if;
      elsif sb(i) /= s(i) then
        return false;
      end if;
      exit when (sb(i) = ' ');
    end loop;
    if s'last < sb'last
     then return false;
     else return true;
    end if;
  end Equal;

  function Classify_Symbol ( sb : Symbol ) return natural32 is
  begin
    if Equal(sb,"time") then
      return 1;
    elsif Equal(sb,"multiplicity") then
      return 2;
    elsif Equal(sb,"err") then
      return 3;
    elsif Equal(sb,"rco") then
      return 4;
    elsif Equal(sb,"res") then
      return 5;
    else
      return 0;
    end if;
  end Classify_Symbol;

  procedure Clear ( L : in out Symbol_List ) is
  begin
     Symbols_Lists.Clear(Symbols_Lists.List(L));
  end Clear;

end Lists_of_Symbols;
