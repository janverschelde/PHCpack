with String_Splitters;
with Characters_and_Numbers;
with Symbol_Table;

package body Sets_of_Unknowns_Strings is

  function to_String ( element : natural32 ) return string is

    standard : constant boolean := (Symbol_Table.Number = 0);

  begin
    if standard then
      declare
        index : constant string := Characters_and_Numbers.nConvert(element);
        res : string(index'first..index'last+1);
      begin
        res(res'first) := 'x';
        res(res'first+1..res'last) := index;
        return res;
      end;
    else
      declare
        sb : Symbol_Table.Symbol;
        res : string(sb'range);
        cnt : integer := sb'first;
      begin
        sb := (sb'range => ' ');
        sb := Symbol_Table.Get(element);
        while sb(cnt) /= ' ' loop
          res(cnt) := sb(cnt);
          cnt := cnt + 1;
        end loop;
        return res(res'first..cnt-1);
      end;
    end if;
  end to_String;

  function to_String ( s : Set ) return string is

    res : String_Splitters.Link_to_String;

  begin
    String_Splitters.Append(res,"{");
    for i in 1..Dimension(s) loop
      if Is_In(s,i) then
        declare
          r : constant string := " " & to_String(i);
        begin
          String_Splitters.Append(res,r);
        end;
      end if;
    end loop;
    String_Splitters.Append(res," }");
    declare
      result : constant string := res.all;
    begin
      String_Splitters.Clear(res);
      return result;
    end;
  end to_String;

  function Parse ( s : string; n : natural32 ) return Set is

    res : Set := Create(n);
    ind : integer;
    standard : constant boolean := (Symbol_Table.Number = 0);

  begin
    ind := s'first;
    while (s(ind) = ' ') loop  -- skip leading spaces
      ind := ind + 1;
      exit when (ind > s'last);
    end loop;
    if ((ind <= s'last) and then (s(ind) = '{')) then 
      -- do not care about trailing spaces
      -- if s(ind) = '{' and s(s'last) = '}' then 
      ind := s'first+1;
      loop
        while s(ind) = ' ' loop
          ind := ind + 1;
        end loop;
        exit when (s(ind) = '}');
        if standard then
          ind := ind + 1; -- skip the x
          declare
            nbr : string(1..integer(n/10)+2);
            cnt : integer := 0;
            index : natural32;
          begin
            while (s(ind) /= ' ') and (s(ind) /= '}') loop
              cnt := cnt + 1;
              nbr(cnt) := s(ind);
              ind := ind + 1;
            end loop;
            index := Characters_and_Numbers.Convert(nbr(1..cnt));
            if index > 0 and index <= n
             then Add(res,index);
            end if;
          end;
        else
          declare
            sb : Symbol_Table.Symbol;
            cnt : integer := 0;
            index : natural32;
          begin
            sb := (sb'range => ' ');
            while (s(ind) /= ' ') and (s(ind) /= '}') loop
              cnt := cnt + 1;
              sb(cnt) := s(ind);
              ind := ind + 1;
            end loop;
            index := Symbol_Table.Get(sb);
            if index > 0 and index <= n
             then Add(res,index);
            end if;
          end;
        end if;
        exit when (s(ind) = '}');
      end loop;
    end if;
    return res;
  end Parse;

end Sets_of_Unknowns_Strings;
