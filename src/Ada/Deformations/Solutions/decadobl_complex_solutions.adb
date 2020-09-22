with unchecked_deallocation;
with DecaDobl_Complex_Numbers_cv;        use DecaDobl_Complex_Numbers_cv;
with Multprec_DecaDobl_Convertors;       use Multprec_DecaDobl_Convertors;
with DecaDobl_Complex_Vectors_cv;        use DecaDobl_Complex_Vectors_cv;

package body DecaDobl_Complex_Solutions is

  use List_of_Solutions;

-- CREATORS :

  function Create ( sl : Solution_List ) return Solution_Array is

    sa : Solution_Array(1..integer32(Length_Of(sl)));

  begin
    if not Is_Null(sl) then
      declare
        i : integer32 := 1;
        temp : Solution_List := sl;
      begin
        while not Is_Null(temp) loop
          sa(i) := new Solution'(Head_Of(temp).all);
          i := i + 1;
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
    return sa;
  end Create;

  function Create ( sa : Solution_Array ) return Solution_List is

    sl : Solution_List;

  begin
    if sa'first <= sa'last then
      declare
        n : constant integer32 := sa(sa'first).n;
        sol : Solution(n) := sa(sa'first).all;
        ls : Link_to_Solution := new Solution'(sol);
        last,tmp : Solution_List;
      begin
        Construct(ls,sl);
        last := sl;
        for i in (sa'first+1)..sa'last loop
          sol := sa(i).all;
          ls := new Solution'(sol);
          Construct(ls,tmp);
          Swap_Tail(last,tmp);
          last := Tail_Of(last);
        end loop;
      end;
    end if;
    return sl;
  end Create;

  function Create ( s : Standard_Complex_Solutions.Solution )
                  return DecaDobl_Complex_Solutions.Solution is

    res : DecaDobl_Complex_Solutions.Solution(s.n);

  begin
    res.t := Standard_to_DecaDobl_Complex(s.t);
    res.m := s.m;
    res.v := Standard_to_DecaDobl_Complex(s.v);
    res.err := Create(s.err);
    res.rco := Create(s.rco);
    res.res := Create(s.res);
    return res;
  end Create;

  function Create ( s : Multprec_Complex_Solutions.Solution )
                  return DecaDobl_Complex_Solutions.Solution is

    res : DecaDobl_Complex_Solutions.Solution(s.n);

  begin
    res.t := Multprec_to_DecaDobl_Complex(s.t);
    res.m := s.m;
    res.v := Multprec_to_DecaDobl_Complex(s.v);
    res.err := to_deca_double(s.err);
    res.rco := to_deca_double(s.rco);
    res.res := to_deca_double(s.res);
    return res;
  end Create;

  function Create ( L : Standard_Complex_Solutions.Solution_List )
                  return DecaDobl_Complex_Solutions.Solution_List is

    res,res_last : DecaDobl_Complex_Solutions.Solution_List;
    tmp : Standard_Complex_Solutions.Solution_List := L;

    use Standard_Complex_Solutions;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Standard_Complex_Solutions.Link_to_Solution
           := Head_Of(tmp);
        ms : constant DecaDobl_Complex_Solutions.Solution(ls.n)
           := Create(ls.all);
      begin
        Append(res,res_last,ms);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( L : Multprec_Complex_Solutions.Solution_List )
                  return DecaDobl_Complex_Solutions.Solution_List is

    res,res_last : DecaDobl_Complex_Solutions.Solution_List;
    tmp : Multprec_Complex_Solutions.Solution_List := L;

    use Multprec_Complex_Solutions;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Multprec_Complex_Solutions.Link_to_Solution
           := Head_Of(tmp);
        ms : constant DecaDobl_Complex_Solutions.Solution(ls.n)
           := Create(ls.all);
      begin
        Append(res,res_last,ms);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

-- COMPARISON and COPYING :

  function Equal ( s1,s2 : Solution; tol : double_float ) return boolean is
  begin
    if (not Equal(s1.t,s2.t)) or else (s1.n /= s2.n) then
      return false;
    else 
      for i in s1.v'range loop
        if i > s2.v'last
         then return false;
        end if;
        if AbsVal(s1.v(i) - s2.v(i)) > tol
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Equal;

  function Equal ( s1,s2 : Solution_List; tol : double_float )
                 return boolean is
  begin
    if Is_Null(s1) and Is_Null(s2) then
      return true;
    elsif Is_Null(s1) or Is_Null(s2) then
      return false;
    else
      declare
        temp1 : Solution_List := s1;
        temp2 : Solution_List := s2;
      begin
        while not Is_Null(temp1) and not Is_Null(s2) loop
          if not Equal(Head_Of(temp1).all,Head_Of(temp2).all,tol) then
            return false;
          else
            temp1 := Tail_Of(temp1);
            temp2 := Tail_Of(temp2);
          end if;
        end loop;
        if Is_Null(temp1) and Is_Null(temp2)
         then return true;
         else return false;
        end if;
      end;
    end if;
  end Equal;

  function Equal ( s1,s2 : Solution_Array; tol : double_float )
                 return boolean is
  begin
    if s1'first /= s2'first then
      return false;
    elsif s1'last /= s2'last then
      return false;
    else
      for i in s1'range loop
        if not Equal(s1(i).all,s2(i).all,tol)
         then return false;
        end if;
      end loop;
    end if;
    return true;
  end Equal;

  procedure Equals ( sols : in out Solution_List; flag : in integer32;
                     tol : in double_float; same : out boolean ) is

    ls1,ls2 : Link_to_Solution;
    tmp1,tmp2 : Solution_List;
    len : constant natural32 := Length_Of(sols);

  begin
    same := false;
    tmp1 := sols;
    for i in 1..len-1 loop
      ls1 := Head_Of(tmp1);
      tmp2 := Tail_Of(tmp1);
      for j in i+1..len loop
        ls2 := Head_Of(tmp2);
        if Equal(ls1.all,ls2.all,tol) then
          same := true;
          ls1.m := flag; Set_Head(tmp1,ls1);
          ls2.m := flag; Set_Head(tmp2,ls2);
        end if;
        tmp2 := Tail_Of(tmp2);
      end loop;
      tmp1 := Tail_Of(tmp1);
    end loop;
  end Equals;

  procedure Equals ( sa : in Solution_Array; x : in Vector; i : in integer32;
                     tol : in double_float; j : in out integer32 ) is

    eq : boolean;

  begin
    while j < i loop
      eq := true;
      for k in x'range loop
        if AbsVal(sa(j).v(k) - x(k)) > tol
         then eq := false;
        end if;
        exit when not eq;
      end loop;
      exit when eq;
      j := j + 1;
    end loop;
  end Equals;

  procedure Copy ( s1 : in Solution; s2 : in out Solution ) is
  begin
    s2.t := s1.t;
    s2.m := s1.m;
    Copy(s1.v,s2.v);
    Copy(s1.err,s2.err);
    Copy(s1.rco,s2.rco);
    Copy(s1.res,s2.res);    
  end Copy;

  procedure Copy ( s1 : in Solution_List; s2 : in out Solution_List ) is
  begin
    Clear(s2);
    if not Is_Null(s1) then
      declare
        temp : Solution_List := s1;
        last : Solution_List;
        n : constant integer32 := Head_Of(s1).n;
        sol : Solution(n) := Head_Of(temp).all;
        ns : Solution(n);
      begin
        Copy(sol,ns);
        declare
          ls : constant Link_to_Solution := new Solution'(ns);
        begin
          Construct(ls,s2);
        end;
        last := s2;
        temp := Tail_Of(temp);
        while not Is_Null(temp) loop
          sol := Head_Of(temp).all;
          declare
            ls : constant Link_to_Solution := new Solution'(sol);
            tmp : Solution_List;
          begin
            Construct(ls,tmp);
            Swap_Tail(last,tmp);
          end;
          last := Tail_Of(last);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Copy;

  procedure Copy ( s1 : in Solution_Array; s2 : in out Solution_Array ) is
  begin
    Clear(s2);
    for i in s1'range loop
      s2(i) := new Solution'(s1(i).all);
    end loop;
  end Copy;

-- SELECTORS :

  function Number ( sols : Solution_List; flag : integer32 )
                  return natural32 is

    res : natural32 := 0;

  begin
    if Is_Null(sols) then
      return res;
    else
      declare
        tmp : Solution_List := sols;
      begin
        while not Is_Null(tmp) loop
          if Head_Of(tmp).m = flag
           then res := res + 1;
          end if;
          tmp := Tail_Of(tmp);
        end loop;
      end;
      return res;
    end if;
  end Number;

  function Is_In ( sols : Solution_List; s : Solution; tol : double_float )
                 return boolean is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      if Equal(Head_Of(tmp).all,s,tol)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( sa : Solution_Array; s : Solution; tol : double_float )
                 return boolean is
  begin
    for i in sa'range loop
      if Equal(sa(i).all,s,tol)
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Retrieve ( sols : Solution_List; p : natural32 )
                    return Link_to_Solution is

    tmp : Solution_List := sols;

  begin
    for i in 1..(p-1) loop
      exit when Is_Null(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    if Is_Null(tmp)
     then return null;
     else return Head_Of(tmp);
    end if;
  end Retrieve;

-- CONSTRUCTORS :

  procedure Push ( s0 : in Solution_List; s1 : in out Solution_List ) is

    tmp : Solution_List := s0;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Construct(ls,s1);
      tmp := Tail_Of(tmp);
    end loop;
  end Push;

  procedure Append ( first,last : in out Solution_List;
                     ls : in Link_to_Solution ) is
  begin
    if Is_Null(first) then
      Construct(ls,first);
      last := first;
    else
      declare
        tmp : Solution_List;
      begin
        Construct(ls,tmp);
        Swap_Tail(last,tmp);
        last := Tail_Of(last);
      end;
    end if;
  end Append;

  procedure Append ( first,last : in out Solution_List; s : in Solution ) is

    ss : Solution(s.n);
    ls : Link_to_Solution;

  begin
    Copy(s,ss);
    ls := new Solution'(ss);
    Append(first,last,ls);
  end Append;

  procedure Add ( sols : in out Solution_List; s : in Solution ) is

    last,temp,tmp : Solution_List;
    ls : constant Link_to_Solution := new Solution'(s);

  begin
    if Is_Null(sols) then
      Construct(ls,sols);
    else
      temp := sols;
      while not Is_Null(temp) loop
        last := temp;
        temp := Tail_Of(temp);
      end loop;
      Construct(ls,tmp);
      Swap_Tail(last,tmp);
    end if;
  end Add;

  procedure Add ( sols : in out Solution_List; s : in Solution;
                  tol : in double_float; other : out natural32 ) is

    last,temp,tmp : Solution_List;
    ls : Link_to_Solution := new Solution'(s);
    s2 : Solution(s.n);
    count : natural32 := 1;

  begin
    other := 0;
    if Is_Null(sols) then
      Construct(ls,sols);
    else
      temp := sols;
      while not Is_Null(temp) loop
        s2 := Head_Of(temp).all;
        if Equal(s,s2,tol) then
          other := count;
          Clear(ls);
          return;
        else
          last := temp;
          temp := Tail_Of(temp);
          count := count + 1;
        end if;
      end loop;
      Construct(ls,tmp);
      Swap_Tail(last,tmp);
    end if;
  end Add;

-- MODIFIERS :

  procedure Change ( sols : in out Solution_List; pos : in natural32;
                     s : in Solution; tol : in double_float;
                     other : out natural32 ) is
  begin
    if pos <= Length_Of(sols) then
      declare
        temp : Solution_List := sols;
        ls : Link_to_Solution;
      begin
        other := 0;
        for i in 1..Length_Of(temp) loop
          ls := Head_Of(temp);
          if i = pos then
            ls.v := s.v;
            ls.m := s.m;
            ls.t := s.t;
            Set_Head(temp,ls);
            return;
          elsif Equal(s,ls.all,tol) then
            other := i;
            return;
          end if;
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Change;

  procedure Set_Continuation_Parameter
               ( sols : in out Solution_List; t : in Complex_Number ) is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Copy(t,ls.t);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Set_Continuation_Parameter;

  procedure Change_Multiplicity
                ( sols : in out Solution_List; pos : in natural32;
                  m : in integer32 ) is
  begin
    if pos <= Length_Of(sols) then
      declare
        temp : Solution_List := sols;
        ls : Link_to_Solution;
      begin
        for i in 1..(pos-1) loop
          temp := Tail_Of(temp);
        end loop;
        ls := Head_Of(temp);
        ls.m := m;
        Set_Head(temp,ls);
      end;
    end if;
  end Change_Multiplicity;

  procedure Remove ( sols : in out Solution_List; pos : in natural32 ) is

    first,second,temp : Solution_List;
    ls : Link_to_Solution;

  begin
    if pos <= Length_Of(sols) then
      if pos = 1 then
        if Is_Null(Tail_Of(sols)) then
          Clear(sols);
        else
          ls := Head_Of(sols);
          Clear(ls);
          sols := Tail_Of(sols);
        end if;
      else
        second := sols;
        for i in 1..(pos-1) loop
          first := second;
          second := Tail_Of(first);
        end loop;
        ls := Head_Of(second);
        Clear(ls);
        temp := Tail_Of(second);
        Swap_Tail(first,temp);
      end if;
    end if;
  end Remove;

  procedure Delete ( sols : in out Solution_List ) is

    continue : boolean;

  begin
    continue := true; -- looking for the first element in sols that stays
    while not Is_Null(sols) and continue loop
      declare
        ls : Link_to_Solution := Head_Of(sols);
      begin
        if To_Be_Removed(ls.m) then
          Clear(ls);
          sols := Tail_Of(sols);
        else
          continue := false;
        end if;
      end;
    end loop;
    if not Is_Null(sols) then -- first element of sols stays in the list
      declare
        first,second : Solution_List;
      begin
        first := sols;
        second := Tail_Of(first);
        while not Is_Null(second) loop
          declare
            ls : Link_to_Solution := Head_Of(second);
            temp : Solution_List;
          begin
            if To_Be_Removed(ls.m) then
              Clear(ls);
              temp := Tail_Of(second);
              Swap_Tail(first,temp);
            end if;
          end;
          first := second;
          second := Tail_Of(first);
        end loop;
      end;
    end if;
  end Delete;
 
  procedure Remove_All ( sols : in out Solution_List; flag : in integer32 ) is

    continue : boolean;

  begin
    continue := true; -- looking for the first element in sols that stays
    while not Is_Null(sols) and continue loop
      declare
        ls : Link_to_Solution := Head_Of(sols);
      begin
        if ls.m = flag then
          Clear(ls);
          sols := Tail_Of(sols);
        else
          continue := false;
        end if;
      end;
    end loop;
    if not Is_Null(sols) then -- first element of s can stay in the list
      declare
        first,second : Solution_List;
      begin
        first := sols;
        second := Tail_Of(first);
        while not Is_Null(second) loop
          declare
            ls : Link_to_Solution := Head_Of(second);
            temp : Solution_List;
          begin
            if ls.m = flag then
              Clear(ls);
              temp := Tail_Of(second);
              Swap_Tail(first,temp);
            end if;
          end;
          first := second;
          second := Tail_Of(first);
        end loop;
      end;
    end if;
  end Remove_All;
    
-- DESTRUCTORS :

  procedure Clear( s : in out Solution ) is
  begin
    Clear(s.err);
    Clear(s.res);
    Clear(s.rco);
    Clear(s.v);
  end Clear;

  procedure Clear ( ls : in out Link_to_Solution ) is

    procedure free is new unchecked_deallocation(Solution,Link_to_Solution);

  begin
    if ls /= null
     then Clear(ls.all);
    end if;
    free(ls);
  end Clear;

  procedure Shallow_Clear ( sl : in out Solution_List ) is
  begin
    List_of_Solutions.Clear(List_of_Solutions.List(sl));
  end Shallow_Clear;

  procedure Deep_Clear ( sl : in out Solution_List ) is

    temp : Solution_List := sl;
    ls : Link_to_Solution;

  begin
    while not Is_Null(temp) loop
      ls := Head_Of(temp);
      Clear(ls);
      temp := Tail_Of(temp);
    end loop;
    Shallow_Clear(sl);
  end Deep_Clear;

  procedure Clear ( sa : in out Solution_Array ) is
  begin
    for i in sa'range loop
      Clear(sa(i));
    end loop;
  end Clear;

  procedure Clear ( s : in out Array_of_Solution_Lists ) is
  begin
    for i in s'range loop
      Deep_Clear(s(i));
    end loop;
  end Clear;

  procedure Clear ( s : in out Link_to_Array_of_Solution_Lists ) is

    procedure free is
      new unchecked_deallocation(Array_of_Solution_Lists,
                                 Link_to_Array_of_Solution_Lists);

  begin
    if s /= null then
      Clear(s.all);
      free(s);
    end if;
  end Clear;

end DecaDobl_Complex_Solutions;
