with unchecked_deallocation;

package body Standard_Continuation_Data is

  use Lists_of_Solu_Info_Arrays;

-- CREATORS :

  function Shallow_Create ( s : Link_to_Solution ) return Solu_Info is

    res : Solu_Info;

  begin
    res.sol := s;
    Init_Info(res);
    res.cora := s.err;
    res.rcond := s.rco;
    res.resa := s.res;
    return res;
  end Shallow_Create;

  function Deep_Create ( s : Solution ) return Solu_Info is

    res : Solu_Info;

  begin
    res.sol := new Solution'(s);
    Init_Info(res);
    res.cora := s.err;
    res.rcond := s.rco;
    res.resa := s.res;
    return res;
  end Deep_Create;

  function Shallow_Create ( s : Solution_Array ) return Solu_Info_Array is

    res : Solu_Info_Array(s'range);

  begin
    for k in res'range loop
      res(k) := Shallow_Create(s(k));
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( s : Solution_Array ) return Solu_Info_Array is

    res : Solu_Info_Array(s'range);

  begin
    for k in res'range loop
      res(k) := Deep_Create(s(k).all);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( s : Solution_List )  return Solu_Info_Array is
  
    res : Solu_Info_Array(1..integer32(Length_Of(s)));
    tmp : Solution_List := s;

  begin
    for k in res'range loop
      res(k) := Shallow_Create(Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Shallow_Create;
   
  function Deep_Create ( s : Solution_List )  return Solu_Info_Array is
   
    res : Solu_Info_Array(1..integer32(Length_Of(s)));
    tmp : Solution_List := s;
 
  begin
    for k in res'range loop
      res(k) := Deep_Create(Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( s : Solu_Info ) return Link_to_Solution is
  begin
    s.sol.err := s.cora;
    s.sol.rco := s.rcond;
    s.sol.res := s.resa;
    return s.sol;
  end Shallow_Create;

  function Deep_Create ( s : Solu_Info ) return Solution is
  begin
    s.sol.err := s.cora;
    s.sol.rco := s.rcond;
    s.sol.res := s.resa;
    return s.sol.all;
  end Deep_Create;

  function Shallow_Create ( s : Solu_Info_Array ) return Solution_Array is

    res : Solution_Array(s'range);

  begin
    for k in s'range loop
      res(k) := Shallow_Create(s(k));
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( s : Solu_Info_Array ) return Solution_Array is

    res : Solution_Array(s'range);

  begin
    for k in s'range loop
      res(k) := new Solution'(Deep_Create(s(k)));
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( s : Solu_Info_Array ) return Solution_List is

    res,res_last : Solution_List;

  begin
    for k in s'range loop
      Append(res,res_last,Shallow_Create(s(k)));
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( s : Solu_Info_Array ) return Solution_List is

    res,res_last : Solution_List;

  begin
    for k in s'range loop
      Append(res,res_last,Deep_Create(s(k)));
    end loop;
    return res;
  end Deep_Create;

  function Create ( s : Solu_Info_Array ) return Solu_Info_Array_List is

    ls : constant Link_to_Solu_Info_Array := new Solu_Info_Array'(s);

  begin
    return Create(ls);
  end Create;

  function Create ( s : Link_to_Solu_Info_Array )
                  return Solu_Info_Array_List is

    res : Solu_Info_Array_List;

  begin
    Construct(s,res);
    return res;
  end Create;

  function Create ( s : Solution_List; size : natural32 )
                  return Solu_Info_Array_List is

    res,res_last : Solu_Info_Array_List;

  begin
    if Length_Of(s) <= size then
      declare
        sa : constant Solu_Info_Array := Deep_Create(s);
      begin
        res := Create(sa);
      end;
    else
      declare
        tmp : Solution_List := s;
        cnt : natural32 := 0;
        buffer : Solu_Info_Array(1..integer32(size));
      begin
        buffer(buffer'first).sol := null;
        while not Is_Null(tmp) loop
          if cnt < size then
            cnt := cnt + 1;
          else
            Append(res,res_last,buffer);
            cnt := 1;
          end if;
          buffer(integer32(cnt)) := Deep_Create(Head_Of(tmp).all);
          tmp := Tail_Of(tmp);
        end loop;
        if cnt > 0
         then Append(res,res_last,buffer(1..integer32(cnt)));
        end if;
      end;
    end if;
    return res;
  end Create;

  procedure Append ( first,last : in out Solu_Info_Array_List;
                     sa : in Solu_Info_Array ) is

    lsa : constant Link_to_Solu_Info_Array := new Solu_Info_Array'(sa);

  begin
    Append(first,last,lsa);
  end Append;

-- CONVERTOR : 

  function Concat ( s : Solu_Info_Array_List ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solu_Info_Array_List := s;

  begin
    while not Is_Null(tmp) loop
      declare
        lsa : constant Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        for i in lsa'range loop
          Append(res,res_last,Shallow_Create(lsa(i)));
        end loop;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Concat;

-- OPERATIONS ON Solu_Info :

  procedure Copy_Info ( s1 : in Solu_Info; s2 : in out Solu_Info ) is
  begin
    s2.corr := s1.corr; s2.cora := s1.cora;
    s2.resr := s1.resr; s2.resa := s1.resa;
    s2.rcond := s1.rcond;  s2.length_path := s1.length_path;
    s2.nstep := s1.nstep; s2.nfail := s1.nfail;
    s2.niter := s1.niter; s2.nsyst := s1.nsyst;
  end Copy_Info;

  procedure Copy_Solu ( s1 : in Solu_Info; s2 : in out Solu_Info ) is
  begin
    Clear(s2.sol);
    s2.sol := new Solution'(s1.sol.all);
  end Copy_Solu;

  procedure Copy ( s1 : in Solu_Info; s2 : in out Solu_Info ) is
  begin
    Copy_Info(s1,s2);
    Copy_Solu(s1,s2);
  end Copy;

  procedure Init_Info ( s : in out Solu_Info ) is
  begin
    s.corr := 0.0; s.cora := 0.0; s.resr := 0.0; s.resa := 0.0; s.rcond := 0.0;
    s.length_path := 0.0;
    s.nstep := 0; s.nfail := 0; s.niter := 0; s.nsyst := 0;
  end Init_Info;

  procedure Add_Info ( s1 : in out Solu_Info; s2 : in Solu_Info ) is
  begin
    s1.nstep := s1.nstep + s2.nstep;
    s1.nfail := s1.nfail + s2.nfail;
    s1.niter := s1.niter + s2.niter;
    s1.nsyst := s1.nsyst + s2.niter;
    s1.length_path := s1.length_path + s2.length_path;
  end Add_Info;

  procedure Update_Info ( s1 : in out Solu_Info; s2 : in Solu_Info ) is
  begin
    s1.corr := s2.corr; s1.cora := s2.cora;
    s1.resr := s2.resr; s1.resa := s2.resa;
    s1.rcond := s2.rcond;
    Add_Info(s1,s2);
  end Update_Info;

-- OPERATIONS ON Solu_Info_Array :

  procedure Copy ( s : in Solu_Info_Array; sa : in out Solution_Array ) is
  begin
    Clear(sa);
    for k in s'range loop
      sa(k) := new Solution'(s(k).sol.all);
    end loop;
  end Copy;

  procedure Copy ( sa : in Solution_Array; s : in out Solu_Info_Array ) is
  begin
    for k in sa'range loop
      Clear(s(k).sol);
      s(k).sol := new Solution'(sa(k).all);
    end loop;
  end Copy;

-- DESTRUCTORS :

  procedure Clear ( s : in out Solu_Info ) is
  begin
    Clear(s.sol);
  end Clear;

  procedure Clear ( s : in out Solu_Info_Array ) is
  begin
    for k in s'range loop
      Clear(s(k));
    end loop;
  end Clear;

  procedure Deep_Clear ( s : in out Link_to_Solu_Info_Array ) is
  begin
    if s /= null
     then Clear(s.all);
          Shallow_Clear(s);
    end if;
  end Deep_Clear;

  procedure Shallow_Clear ( s : in out Link_to_Solu_Info_Array ) is

    procedure free is
      new unchecked_deallocation(Solu_Info_Array,Link_to_Solu_Info_Array);

  begin
    free(s);
  end Shallow_Clear;

  procedure Deep_Clear ( s : in out Solu_Info_Array_List ) is

    tmp : Solu_Info_Array_List := s;

  begin
    while not Is_Null(tmp) loop
      declare
        lsa : Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        Deep_Clear(lsa);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Solu_Info_Arrays.Clear(Lists_of_Solu_Info_Arrays.List(s));
  end Deep_Clear;

  procedure Shallow_Clear ( s : in out Solu_Info_Array_List ) is
  begin
    Lists_of_Solu_Info_Arrays.Clear(Lists_of_Solu_Info_Arrays.List(s));
  end Shallow_Clear;

end Standard_Continuation_Data;
