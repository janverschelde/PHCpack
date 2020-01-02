with unchecked_deallocation;
with Standard_Integer_Numbers;             use Standard_Integer_Numbers;

package body Standard_Solution_Array_Lists is

-- CREATORS :

  function Create ( s : Solution_Array ) return Solution_Array_List is

    ls : constant Link_to_Solution_Array := new Solution_Array'(s);

  begin
    return Create(ls);
  end Create;

  function Create ( s : Link_to_Solution_Array )
                  return Solution_Array_List is

    res : Solution_Array_List;

  begin
    Construct(s,res);
    return res;
  end Create;

  function Create ( s : Solution_List; size : natural32 )
                  return Solution_Array_List is

    res,res_last : Solution_Array_List;

  begin
    if Length_Of(s) <= size then
      declare
        sa : constant Solution_Array := Create(s);
      begin
        res := Create(sa);
      end;
    else
      declare
        tmp : Solution_List := s;
        cnt : natural32 := 0;
        buffer : Solution_Array(1..integer32(size));
      begin
        while not Is_Null(tmp) loop
          if cnt < size
           then cnt := cnt + 1;
           else Append(res,res_last,buffer);
                cnt := 1;
          end if;
          buffer(integer32(cnt)) := Head_Of(tmp);
          tmp := Tail_Of(tmp);
        end loop;
        if cnt > 0
         then Append(res,res_last,buffer(1..integer32(cnt)));
        end if;
      end;
    end if;
    return res;
  end Create;

  procedure Append ( first,last : in out Solution_Array_List;
                     sa : in Solution_Array ) is

    lsa : constant Link_to_Solution_Array := new Solution_Array'(sa);

  begin
    Append(first,last,lsa);
  end Append;

-- CONVERTOR : 

  function Concat ( s : Solution_Array_List ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_Array_List := s;

  begin
    while not Is_Null(tmp) loop
      declare
        lsa : constant Link_to_Solution_Array := Head_Of(tmp);
      begin
        for i in lsa'range loop
          Append(res,res_last,lsa(i));
        end loop;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Concat;

-- DESTRUCTORS :

  procedure Deep_Clear ( s : in out Link_to_Solution_Array ) is
  begin
    if s /= null
     then Clear(s.all);
          Shallow_Clear(s);
    end if;
  end Deep_Clear;

  procedure Shallow_Clear ( s : in out Link_to_Solution_Array ) is

    procedure free is
      new unchecked_deallocation(Solution_Array,Link_to_Solution_Array);

  begin
    free(s);
  end Shallow_Clear;

  procedure Deep_Clear ( s : in out Solution_Array_List ) is

    tmp : Solution_Array_List := s;

  begin
    while not Is_Null(tmp) loop
      declare
        lsa : Link_to_Solution_Array := Head_Of(tmp);
      begin
        Deep_Clear(lsa);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Solution_Arrays.Clear(Lists_of_Solution_Arrays.List(s));
  end Deep_Clear;

  procedure Shallow_Clear ( s : in out Solution_Array_List ) is
  begin
    Lists_of_Solution_Arrays.Clear(Lists_of_Solution_Arrays.List(s));
  end Shallow_Clear;

end Standard_Solution_Array_Lists;
