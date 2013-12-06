with unchecked_deallocation;

package body Trees_of_Vectors is

-- SELECTORS :

  function Is_In ( tv : Tree_of_Vectors; v : Vector ) return boolean is

    tmp : Tree_of_Vectors;
    d2 : Link_to_Vector;

  begin
    tmp := tv;
    while not Is_Null(tmp) loop
      d2 := Head_Of(tmp).d;
      if Equal(d2.all,v)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( tv : Tree_of_Vectors; v : Link_to_Vector ) return boolean is
  begin
    if v /= null
     then return Is_In(tv,v.all);
     else return false;
    end if;
  end Is_In;

  procedure Iterator ( tv : in Tree_of_Vectors ) is

    tmp : Tree_of_Vectors;
    cont : boolean;

  begin
    tmp := tv;
    while not Is_Null(tmp) loop
      Process(Head_Of(tmp),cont);
      exit when not cont;
      tmp := Tail_Of(tmp);
    end loop;
  end Iterator;

-- DESTRUCTORS :

  procedure Clear ( nd : in out node ) is
  begin
    Clear(nd.d);
    Clear(nd.ltv);
  end Clear;

  procedure Clear ( tv : in out Tree_of_Vectors ) is

    tmp : Tree_of_Vectors;

  begin
    tmp := tv;
    while not Is_Null(tmp) loop
      declare
        nd : node := Head_Of(tmp);
      begin
        Clear(nd);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Link_to_Vector_Trees.Clear(Link_to_Vector_Trees.List(tv));
  end Clear;

  procedure Clear ( ltv : in out Link_to_Tree_of_Vectors ) is

    procedure free is 
      new unchecked_deallocation(Tree_of_Vectors,Link_to_Tree_of_Vectors);

  begin
    if not (ltv = null)
     then Clear(ltv.all); free(ltv);
    end if;
  end Clear;

end Trees_of_Vectors;
