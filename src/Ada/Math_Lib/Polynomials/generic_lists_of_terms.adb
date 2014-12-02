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

  procedure Append ( first,last : in out Term_List; t : in Term ) is
  begin
    List_of_Terms.Append(List_of_Terms.List(first),List_of_Terms.List(last),t);
  end Append;

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

-- COPYING :

  procedure Copy ( p : in Term_List; q : in out Term_List ) is

    tmp : Term_List := p;
    q_last : Term_List;
    t : Term;

  begin
    Clear(q); q_last := q;
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
    List_of_Terms.Clear(List_of_Terms.List(p));
  end Shallow_Clear;

  procedure Shallow_Clear ( p : in out Array_of_Term_Lists ) is
  begin
    for i in p'range loop
      Shallow_Clear(p(i));
    end loop;
  end Shallow_Clear;

end Generic_Lists_of_Terms;
