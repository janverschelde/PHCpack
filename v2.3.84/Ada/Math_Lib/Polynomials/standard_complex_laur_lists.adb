with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Standard_Complex_Laur_Lists is

-- CONSTRUCTORS :

  function Create ( p : Poly ) return Poly_List is

    res : Poly_List;

  begin
    Construct(p,res);
    return res;
  end Create;

  function Create_Copy ( p : Poly ) return Poly_List is

    res : Poly_List;
    q : Poly;

  begin
    Standard_Complex_Laurentials.Copy(p,q);
    Construct(q,res);
    return res;
  end Create_Copy;

  function Create ( p : Laur_Sys ) return Poly_List is

    res,res_last : Poly_List;

  begin
    for i in p'range loop
      Append(res,res_last,p(i));
    end loop;
    return res;
  end Create;

  function Create_Copy ( p : Laur_Sys ) return Poly_List is

    res,res_last : Poly_List;

  begin
    for i in p'range loop
      Append_Copy(res,res_last,p(i));
    end loop;
    return res;
  end Create_Copy;

  function Create ( p : Poly_List ) return Laur_Sys is

    res : Laur_Sys(1..integer32(Length_Of(p)));
    tmp : Poly_List := p;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create_Copy ( p : Poly_List ) return Laur_Sys is

    res : Laur_Sys(1..integer32(Length_Of(p)));
    tmp : Poly_List := p;

  begin
    for i in res'range loop
      Standard_Complex_Laurentials.Copy(Head_Of(tmp),res(i));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create_Copy;

  procedure Append ( first,last : in out Poly_List; p : in Poly ) is
  begin
    if Is_Null(first) then 
      Construct(p,first);
      last := first;
    else
      declare
        tmp : Poly_List;
      begin
        Construct(p,tmp);
        Swap_Tail(last,tmp);
        last := Tail_Of(last);
      end;
    end if;
  end Append;

  procedure Append_Copy ( first,last : in out Poly_List; p : in Poly ) is

    q : Poly;

  begin
    Standard_Complex_Laurentials.Copy(p,q);
    Append(first,last,q);
  end Append_Copy;

  procedure Concatenate ( first,last : in out Poly_List;
                          pl : in Poly_List ) is

    tmp : Poly_List := pl;

  begin
    while not Is_Null(tmp) loop
      Append(first,last,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end Concatenate;

-- DESTRUCTORS :

  procedure Shallow_Clear ( pl : in out Poly_List ) is
  begin
    List_of_Laurentials.Clear(List_of_Laurentials.List(pl));
  end Shallow_Clear;

  procedure Deep_Clear ( pl : in out Poly_List ) is

    tmp : Poly_List := pl;
    p : Poly;

  begin
    while not Is_Null(tmp) loop
      p := Head_Of(tmp);
      Standard_Complex_Laurentials.Clear(p);
      tmp := Tail_Of(tmp);
    end loop;
    List_of_Laurentials.Clear(List_of_Laurentials.List(pl));
  end Deep_Clear;

end Standard_Complex_Laur_Lists;
