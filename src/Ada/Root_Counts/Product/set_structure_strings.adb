with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with String_Splitters;
with Symbol_Table;
with Set_Structure;

package body Set_Structure_Strings is

-- WRITE STRING REPRESENTATIONS :

  function to_String ( i,j : in natural32 ) return string is
  begin
    if Set_Structure.Empty then
      return "";
    else
      declare
        res : String_Splitters.Link_to_String;
        n : constant natural32 := Set_Structure.Dimension;
        sb : Symbol_Table.Symbol;
        ind : integer;
      begin
        String_Splitters.Append(res,"{ ");
        for k in 1..n loop
          if Set_Structure.Is_In(i,j,k) then
            sb := Symbol_Table.Get(k);
            ind := sb'first;
            while sb(ind) /= ' ' loop
              ind := ind + 1;
            end loop;
            String_Splitters.Append(res,sb(sb'first..ind));
          end if;
        end loop;
        String_Splitters.Append(res,"}");
        declare
          result : constant string := res.all;
        begin
          String_Splitters.Clear(res);
          return result;
        end;
      end;
    end if;
  end to_String;

  function to_String ( i : in natural32 ) return string is

    res : String_Splitters.Link_to_String;

  begin
    if Set_Structure.Empty then
      return "";
    else
      for j in 1..Set_Structure.Number_of_Sets(i) loop
        declare
          jset : constant string := to_String(i,j);
        begin
          String_Splitters.Append(res,jset);
        end;
      end loop;
      declare
        result : constant string := res.all;
      begin
        String_Splitters.Clear(res);
        return result;
      end;
    end if;
  end to_String;

  function to_String return string is

    res : String_Splitters.Link_to_String;

  begin
    if Set_Structure.Empty then
      return "";
    else
      for i in 1..Set_Structure.Dimension loop
        declare
          equ : constant string := to_String(i);
        begin
          String_Splitters.Append(res,equ & ';');
        end;
      end loop;
      declare
        result : constant string := res.all;
      begin
        String_Splitters.Clear(res);
        return result;
      end;
    end if;
  end to_String;

-- PARSE STRING REPRESENTATIONS :

  function Number_of_Semicolons ( s : string ) return natural32 is

  -- DESCRIPTION :
  --   Counts the number of semicolons in the string.

    res : natural32 := 0;

  begin
    for i in s'range loop
      if s(i) = ';'
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Semicolons;

  function Number_of_Sets
             ( s : string ) return Standard_Natural_Vectors.Vector is

    n : constant natural32 := Number_of_Semicolons(s);
    res : Standard_Natural_Vectors.Vector(1..integer32(n));
    ind : integer32 := res'first;

  begin
    res := (res'range => 0);
    for i in s'range loop
      if s(i) = ';' then
        ind := ind + 1;
      elsif s(i) = '{' then
        res(ind) := res(ind) + 1;
      end if;
    end loop;
    return res;
  end Number_of_Sets;

  procedure Parse ( s : in string ) is

    ns : constant Standard_Natural_Vectors.Vector := Number_of_Sets(s);
    startind,nextind : integer := s'first;

  begin
    Set_Structure.Clear;
    Set_Structure.Init(ns);
    for i in 1..natural32(ns'last) loop
      while s(nextind) /= ';' loop
        nextind := nextind + 1;
        exit when (nextind > s'last);
      end loop;
      Parse(s(startind..nextind-1),i);
      nextind := nextind + 1;
      exit when (nextind > s'last);
      startind := nextind;
    end loop;
  end Parse;

  procedure Parse ( s : in string; i : in natural32 ) is

    j,k : natural32 := 0;
    sb : Symbol_Table.Symbol;
    ind : integer := s'first;
    sbind : integer;

  begin
    while ind <= s'last loop
      if s(ind) = ' ' then       -- skip space
        ind := ind + 1;
      elsif s(ind) = '{' then    -- start of new set
        j := j + 1;
        ind := ind + 1;
      elsif s(ind) = '}' then    -- ignore end of set
        ind := ind + 1;
      else                       -- parse symbol
        sb := (sb'range => ' ');
        sb(sb'first) := s(ind);
        ind := ind + 1;
        sbind := sb'first+1;
        while s(ind) /= ' ' and s(ind) /= '}' loop
          sb(sbind) := s(ind);
          ind := ind + 1;
          sbind := sbind + 1;
          exit when (ind > s'last);
        end loop;
        k := Symbol_Table.get(sb);
        if k /= 0
         then Set_Structure.Add(i,j,k);
        end if;
      end if;
    end loop;
  end Parse;

end Set_Structure_Strings;
