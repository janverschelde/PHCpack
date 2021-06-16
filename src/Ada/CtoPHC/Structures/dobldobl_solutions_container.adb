package body DoblDobl_Solutions_Container is

-- INTERNAL DATA :

  first,last,current : Solution_List;
  cursor : natural32 := 0; -- index to the current solution

-- OPERATIONS :

  procedure Initialize ( sols : in Solution_List ) is

    tmp : Solution_List := sols;
 
  begin
    for i in 1..Length_Of(sols) loop
      Append(first,last,Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    current := first;
    cursor := 1;
  end Initialize;

  function Length return natural32 is
  begin
    return Length_Of(first);
  end Length;

  function Dimension return natural32 is
  begin
    if Is_Null(first)
     then return 0;
     else return natural32(Head_Of(first).n);
    end if;
  end Dimension;

  function Retrieve return Solution_List is
  begin
    return first;
  end Retrieve;

  procedure Retrieve ( k : in natural32; s : out Solution;
                       fail : out boolean ) is

    ls : Link_to_Solution;

  begin
    Retrieve(k,ls,fail);
    if not fail
     then s := ls.all;
    end if;
  end Retrieve;

  procedure Retrieve ( k : in natural32; s : out Link_to_Solution;
                       fail : out boolean ) is

    tmp : Solution_List := first;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k then
        fail := false;
        s := Head_Of(tmp);
        return;
      else
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    fail := true;
  end Retrieve;

  procedure Retrieve_Next_Initialize is
  begin
    current := first;
    if Is_Null(current)
     then cursor := 0; -- empty solution list
     else cursor := 1; -- index to the first solution
    end if;
  end Retrieve_Next_Initialize;

  procedure Retrieve_Next ( s : out Link_to_Solution; k : out natural32 ) is
  begin
    if Is_Null(current) then
      k := 0;
    else
      s := Head_Of(current);
      k := cursor;
      current := Tail_Of(current);
      cursor := cursor + 1;
    end if;
  end Retrieve_Next;

  procedure Retrieve_Current ( s : out Link_to_Solution; k : out natural32 ) is
  begin
    if Is_Null(current) then
      k := 0;
    else
      s := Head_Of(current);
      k := cursor;
    end if;
  end Retrieve_Current;

  procedure Move_Current ( k : out natural32 ) is
  begin
    if Is_Null(current) then
      k := 0;
    else
      current := Tail_Of(current);
      if Is_Null(current)
       then cursor := 0;
       else cursor := cursor + 1;
      end if;
      k := cursor;
    end if;
  end Move_Current;

  procedure Replace ( k : in natural32; s : in Solution;
                      fail : out boolean ) is
	  
    tmp : Solution_List := first;
    ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k then
        fail := false;
        ls := Head_Of(tmp);
        ls.t := s.t;
        ls.m := s.m;
        ls.v := s.v;
        ls.err := s.err;
        ls.rco := s.rco;
        ls.res := s.res;
        Set_Head(tmp,ls);
        return;
      else
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    fail := true;
  end Replace;

  procedure Replace ( k : in natural32; s : in Link_to_Solution;
                      fail : out boolean ) is
	  
    tmp : Solution_List := first;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k then
        fail := false;
        Set_Head(tmp,s);
        return;
      else
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    fail := true;
  end Replace;

  procedure Append ( s : in Solution ) is
  begin
    if Is_Null(first) then
      Append(first,last,s);
      current := first;
      cursor := 1;
    else
      Append(first,last,s);
    end if;
  end Append;

  procedure Append ( s : in Link_to_Solution ) is
  begin
    if Is_Null(first) then
      Append(first,last,s);
      current := first;
      cursor := 1;
    else
      Append(first,last,s);
    end if;
  end Append;

  procedure Clear is
  begin
    Clear(first);
    last := first;
    current := first;
    cursor := 0;
  end Clear;

end DoblDobl_Solutions_Container;
