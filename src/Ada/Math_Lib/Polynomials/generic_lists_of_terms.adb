with Standard_Natural_Vectors;

package body Generic_Lists_of_Terms is 

-- CONSTRUCTORS :

  function Create ( p : Poly ) return Term_List is

    res,res_last : Term_List;

    procedure Append_Term ( t : in Term; continue : out boolean ) is
    begin
      Append(res,res_last,t);
      continue := True;
    end Append_Term;
    procedure Append_Terms is new Visiting_Iterator(Append_Term);

  begin
    Append_Terms(p);
    return res;
  end Create;

  function Create ( t : Term_List ) return Poly is

    res : Poly := Null_Poly;
    tmp : Term_List := t;

  begin
    while not Is_Null(tmp) loop
      Add(res,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  procedure Append ( first,last : in out Term_List; t : in Term ) is

    tt : Term;

  begin
    Copy(t.cf,tt.cf);
    Copy(t.dg,tt.dg);
    if Is_Null(first) then
      Construct(tt,first);
      last := first;
    else
      declare
        tmp : Term_List;
      begin
        Construct(tt,tmp);
        Swap_Tail(last,tmp);
        last := Tail_Of(last);
      end;
    end if;
  end Append;

  procedure Merge_Append ( first,last : in out Term_List; t : in Term ) is

    tt : Term;

  begin
    if Is_Null(first) then
      Copy(t.cf,tt.cf);
      Copy(t.dg,tt.dg);
      Construct(tt,first);
      last := first;
    else
      declare
        tmp : Term_List := first;
        ttt : Term;
        newtmp : Term_List;
      begin
        while not Is_Null(tmp) loop  -- search for matching exponents
          ttt := Head_Of(tmp);
          if Standard_Natural_Vectors.Equal(ttt.dg.all,t.dg.all) then
            Add(ttt.cf,t.cf);
            Set_Head(tmp,ttt); return;
          else
            tmp := Tail_Of(tmp);
          end if;
        end loop;
        if Is_Null(tmp) then
          Copy(t.cf,tt.cf);
          Copy(t.dg,tt.dg);
          Construct(tt,newtmp);
          Swap_Tail(last,newtmp);
          last := Tail_Of(last);
        end if;
      end;
    end if;
  end Merge_Append;

  procedure Concat ( first,last : in out Term_List; t : in Term_List ) is

    tmp : Term_List := t;
    tt : Term;

  begin
    while not Is_Null(tmp) loop
      tt := Head_Of(tmp);
      Append(first,last,tt);
      tmp := Tail_Of(tmp);
    end loop;
  end Concat;

  procedure Merge_Concat ( first,last : in out Term_List;
                           t : in Term_List ) is

    tmp : Term_List := t;
    tt : Term;

  begin
    while not Is_Null(tmp) loop
      tt := Head_Of(tmp);
      Merge_Append(first,last,tt);
      tmp := Tail_Of(tmp);
    end loop;
  end Merge_Concat;

-- COPYING :

  procedure Copy ( p : in Term_List; q,q_last : in out Term_List ) is

    tmp : Term_List := p;
    t : Term;

  begin
    if not Is_Null(q)
     then Clear(q); q_last := q;
    end if;
    while not Is_Null(tmp) loop
      t := Head_Of(tmp);
      Append(q,q_last,t);
      tmp := Tail_Of(tmp);
    end loop;
  end Copy;

-- SELECTORS :

  function Is_In ( p : in Term_List; t : in Term ) return boolean is

    tmp : Term_List := p;
    tp : Term;

  begin
    while not Is_Null(tmp) loop
      tp := Head_Of(tmp);
      if Equal(tp,t)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  procedure Iterator ( p : in Term_List ) is

    tmp : Term_List := p;
    t : Term;
    cont : boolean := true;

  begin
    while not Is_Null(tmp) loop
      t := Head_Of(tmp);
      process(t,cont);
      exit when not cont;
      tmp := Tail_Of(tmp);
    end loop;
  end Iterator;

-- DESTRUCTORS :

  procedure Clear ( p : in out Term_List ) is

    tmp : Term_List := p;
    t : Term;

  begin
    while not Is_Null(tmp) loop
      t := Head_Of(tmp);
      Clear(t);
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(p);
  end Clear;

  procedure Clear ( p : in out Array_of_Term_Lists ) is
  begin
    for i in p'range loop
      Clear(p(i));
    end loop;
  end Clear;

  procedure Shallow_Clear ( p : in out Term_List ) is
  begin
    if not Is_Null(p)
     then List_of_Terms.Clear(List_of_Terms.List(p));
    end if;
  end Shallow_Clear;

  procedure Shallow_Clear ( p : in out Array_of_Term_Lists ) is
  begin
    for i in p'range loop
      Shallow_Clear(p(i));
    end loop;
  end Shallow_Clear;

end Generic_Lists_of_Terms;
