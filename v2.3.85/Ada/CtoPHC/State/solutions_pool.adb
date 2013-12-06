package body Solutions_Pool is

-- INTERNAL DATA :

  size_pool : integer32;
  pool_first,pool_last : Link_to_Array_of_Solution_Lists;

-- CREATORS :

  procedure Initialize ( n : in integer32 ) is
  begin
    size_pool := n;
    pool_first := new Array_of_Solution_Lists(1..n);
    pool_last := new Array_of_Solution_Lists(1..n);
  end Initialize;

  procedure Initialize ( k : in integer32; sols : in Solution_List ) is
 
    tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    if k > 0 and k <= size_pool then
      tmp := sols;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        Append(pool_first(k),pool_last(k),ls.all);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Initialize;

-- SELECTORS :

  function Size return natural32 is
  begin
    return natural32(size_pool);
  end Size;

  function Length ( k : in integer32 ) return natural32 is
  begin
    if k > size_pool or k < 1
     then return 0;
     else return Length_Of(pool_first(k));
    end if;
  end Length;

  function Dimension ( k : in integer32 ) return natural32 is
  begin
    if k > size_pool or k < 1 then
      return 0;
    elsif Is_Null(pool_first(k)) then
      return 0;
    else
      return natural32(Head_Of(pool_first(k)).n);
    end if;
  end Dimension;

  function Retrieve ( k : in integer32 ) return Solution_List is

    res : Solution_List;

  begin
    if k > 0 and k <= size_pool
     then res := pool_first(k);
    end if;
    return res;
  end Retrieve;

  procedure Retrieve ( k : in integer32; i : in natural32;
                       s : out Solution; fail : out boolean ) is

    tmp : Solution_List;
    cnt : natural32 := 0;
    ls : Link_to_Solution;

  begin
    if k < 1 or k > size_pool then
      fail := true;
    elsif i > Length_Of(pool_first(k)) then
      fail := true;
    else
      tmp := pool_first(k);
      while not Is_Null(tmp) loop
        cnt := cnt + 1;
        if cnt = i then
          fail := false;
          ls := Head_Of(tmp);
          s := ls.all;
          return;
        else
          tmp := Tail_Of(tmp);
        end if;
      end loop;
    end if;
    fail := true;
  end Retrieve;

  procedure Retrieve ( k : in integer32; i : in natural32;
                       s : out Link_to_Solution; fail : out boolean ) is

    tmp : Solution_List;
    cnt : natural32 := 0;

  begin
    if k < 1 or k > size_pool then
      fail := true;
    elsif i > Length_Of(pool_first(k)) then
      fail := true;
    else
      tmp := pool_first(k);
      while not Is_Null(tmp) loop
        cnt := cnt + 1;
        if cnt = i then
          fail := false;
          s := Head_Of(tmp);
          return;
        else
          tmp := Tail_Of(tmp);
        end if;
      end loop;
    end if;
    fail := true;
  end Retrieve;

  procedure Replace ( k : in integer32; i : in natural32;
                      s : in Solution; fail : out boolean ) is

    tmp : Solution_List;
    cnt : natural32 := 0;
    ls : Link_to_Solution;

  begin
    if k < 1 or k > size_pool then
      fail := true;
    elsif i > Length_Of(pool_first(k)) then
      fail := true;
    else
      tmp := pool_first(k);
      while not Is_Null(tmp) loop
        cnt := cnt + 1;
        if cnt = i then
          fail := false;
          ls := Head_Of(tmp);
          ls.t := s.t; ls.m := s.m; ls.v := s.v;
          ls.err := s.err; ls.rco := s.rco; ls.res := s.res;
          Set_Head(tmp,ls);
          return;
        else
          tmp := Tail_Of(tmp);
        end if;
      end loop;
    end if;
    fail := true;
  end Replace;

  procedure Replace ( k : in integer32; i : in natural32;
                      s : in Link_to_Solution; fail : out boolean ) is

    tmp : Solution_List;
    cnt : natural32 := 0;

  begin
    if k < 1 or k > size_pool then
      fail := true;
    elsif i > Length_Of(pool_first(k)) then
      fail := true;
    else
      tmp := pool_first(k);
      while not Is_Null(tmp) loop
        cnt := cnt + 1;
        if cnt = i
         then fail := false; Set_Head(tmp,s); return;
         else tmp := Tail_Of(tmp);
        end if;
      end loop;
    end if;
    fail := true;
  end Replace;

  procedure Append ( k : in integer32; s : in Solution ) is
  begin
    if k > 0 and k <= size_pool
     then Append(pool_first(k),pool_last(k),s);
    end if;
  end Append;

  procedure Append ( k : in integer32; s : in Link_to_Solution ) is
  begin
    if k > 0 and k <= size_pool
     then Append(pool_first(k),pool_last(k),s);
    end if;
  end Append;

-- DESTRUCTORS :

  procedure Clear ( k : in integer32 ) is
  begin
    if k > 0 and k <= size_pool then
      Clear(pool_first(k));
      pool_last(k) := pool_first(k);
    end if;
  end Clear;

  procedure Clear is
  begin
    size_pool := 0;
  end Clear;

begin
  size_pool := 0;
end Solutions_Pool;
