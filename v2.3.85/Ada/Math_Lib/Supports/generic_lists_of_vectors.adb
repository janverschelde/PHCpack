with Standard_Integer_Numbers;            use Standard_Integer_Numbers;

package body Generic_Lists_of_Vectors is

-- CONSTRUCTORS :

  function Deep_Create ( v : VecVec ) return List is

    res,res_last : List;

  begin
    for i in v'range loop
      Append(res,res_last,v(i).all);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( v : VecVec ) return List is

    res,res_last : List;

  begin
    for i in v'range loop
      Append(res,res_last,v(i));
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( L : List ) return VecVec is

    res : VecVec(1..integer32(Length_Of(L)));
    tmp : List := L;

  begin
    for i in res'range loop
      declare
        v : constant Vectors.Vector := Head_Of(tmp).all;
      begin
        res(i) := new vector'(v);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( L : List ) return VecVec is

    res : VecVec(1..integer32(Length_Of(L)));
    tmp : List := L;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Shallow_Create;

  procedure Copy ( L1 : in List; L2 : in out List ) is

    tmp,L2_last : List;
    lv : Link_to_Vector;

  begin
    Deep_Clear(L2);
    tmp := L1;
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      Append(L2,L2_last,lv.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Copy;

  procedure Append ( first,last : in out List; v : in Vector ) is

    lv : Link_to_Vector := new Vector'(v);

  begin
    if Is_Null(first) then
      Construct(lv,first);
      last := first;
    else
      declare
        tmp : List;
      begin
        Construct(lv,tmp);
        Swap_Tail(last,tmp);
        last := Tail_Of(last);
      end;
    end if;
  end Append;

  procedure Append_Diff ( first,last : in out List; v : in Vector ) is
  begin
    if not Is_In(first,v)
     then Append(first,last,v);
    end if;
  end Append_Diff;

  procedure Append_Diff ( first,last : in out List; v : in Link_to_Vector ) is
  begin
    if v /= null and then not Is_In(first,v)
     then Append(first,last,v);
    end if;
  end Append_Diff;

  procedure Deep_Concat ( first,last : in out List; L : in List ) is

    tmp : List;
    lv : Link_to_Vector;

  begin
    if not Is_Null(L) then
      tmp := L;
      while not Is_Null(tmp) loop
        lv := Head_Of(tmp);
        Append(first,last,lv.all);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Deep_Concat;

  procedure Shallow_Concat ( first,last : in out List; L : in List ) is
  begin
    Concat(first,last,L);
  end Shallow_Concat;

  procedure Deep_Concat_Diff ( first,last : in out List; L : in List ) is

    tmp : List;
    lv : Link_to_Vector;

  begin
    if not Is_Null(L) then
      tmp := L;
      while not Is_Null(tmp) loop
        lv := Head_Of(tmp);
        Append_Diff(first,last,lv.all);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Deep_Concat_Diff;

  procedure Shallow_Concat_Diff ( first,last : in out List; L : in List ) is

    tmp : List;
    lv : Link_to_Vector;

  begin
    if not Is_Null(L) then
      tmp := L;
      while not Is_Null(tmp) loop
        lv := Head_Of(tmp);
        Append_Diff(first,last,lv);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Shallow_Concat_Diff;

  procedure Remove ( L : in out List; x : in Vector ) is

    lpt : Link_to_Vector;
    found : boolean;
    L1,L2 : List;

  begin
    if not Is_Null(L) then
      lpt := Head_Of(L);
      if lpt.all = x then
        Clear(lpt);
        L := Tail_Of(l);
      else 
        found := false;
        L1 := L;
        L2 := Tail_Of(L1);
        while not Is_Null(L2) loop
          lpt := Head_Of(L2);
          found := (lpt.all = x);
          exit when found;
          L1 := L2;
          L2 := Tail_Of(L1);
        end loop;
        if found then
          Clear(lpt);
          L2 := Tail_Of(L2);
          Swap_Tail(L1,L2);
        end if;
       end if;
    end if;
  end Remove;

  procedure Remove ( L : in out List; x : in Link_to_Vector ) is
  begin
    if x /= null
     then Remove(L,x.all);
    end if;
  end Remove;

  procedure Swap_to_Front ( L : in out List; x : in Vector ) is

    first : Link_to_Vector;
    pt : Link_to_Vector;
    tmp : List;
    done : boolean := false;

  begin
    if not Is_Null(L) then
      first := Head_Of(L);
      if first.all /= x then
        tmp := Tail_Of(L);
        while not Is_Null(tmp) loop
          pt := Head_Of(tmp);
          if pt.all = x then
            Set_Head(tmp,first);
            Set_Head(l,pt);
            done := true;
          end if;
          exit when done;
          tmp := Tail_Of(tmp);
        end loop;
      end if;
    end if;
  end Swap_to_Front;

  procedure Swap_to_Front ( L : in out List; x : in Link_to_Vector ) is
  begin
    if x /= null
     then Swap_to_Front(L,x.all);
    end if;
  end Swap_to_Front;

-- SELECTORS :

  function Is_In ( L : List; v : Vector ) return boolean is

    tmp : List;
    v2 : Link_to_Vector;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      v2 := Head_Of(tmp);
      if Equal(v2.all,v)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( L : List; v : Link_to_Vector ) return boolean is
  begin
    if v = null
     then return false;
     else return Is_In(L,v.all);
    end if;
  end Is_In;

  function Sub_List ( L1,L2 : List ) return boolean is

    tmp : List := L1;

  begin
    while not Is_Null(tmp) loop
      if not Is_In(L2,Head_Of(tmp))
       then return false;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return true;
  end Sub_List;

  function Equal ( L1,L2 : List ) return boolean is
  begin
    if not Sub_List(L1,L2) then
      return false;
    elsif not Sub_List(L2,L1) then
      return false;
    else
      return true;
    end if;
  end Equal;

-- DESTRUCTORS :

  procedure Deep_Clear ( L : in out List ) is

    tmp : List;
    v : Link_to_Vector;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      v := Head_Of(tmp);
      Clear(v);
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(l);
  end Deep_Clear;

  procedure Shallow_Clear ( L : in out List ) is
  begin
    Vector_Lists.Clear(Vector_Lists.List(L));
  end Shallow_Clear;

end Generic_Lists_of_Vectors;
