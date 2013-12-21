with Integer32_Vectors_Utilities;        use Integer32_Vectors_Utilities;

package body Transforming_Integer32_Vector_Lists is

  procedure Shift ( L : in out List; v : in Vector ) is

    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      declare
        lv : constant Link_to_Vector := Head_Of(tmp);
      begin
        lv.all := lv.all - v;
        Set_Head(tmp,lv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Shift;

  procedure Shift ( L : in out List; v : in Link_to_Vector ) is
  begin
    if v /= null
     then Shift(L,v.all);
    end if;
  end Shift;

  function Shift ( L : List; v : Vector ) return List is

    tmp,res,res_last : List;
    v1 : Vector(v'range);

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      v1 := Head_Of(tmp).all;
      declare
        v2 : constant Vector(v1'range) := v1 - v;
      begin
        Append(res,res_last,v2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Shift;

  function Shift ( L : List; v : Link_to_Vector ) return List is
  begin
    if v = null then
      declare
        res : List;
      begin
        Copy(L,res);
        return res;
      end;
    else
      return Shift(L,v.all);
    end if;
  end Shift;

  function "*"( L : List; t : Transfo ) return List is
  begin
    return t*L;
  end "*";

  function "*"( t : Transfo; L : List ) return List is

    tmp,res,res_last : List;
    d,td : Link_to_Vector;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      d := Head_Of(tmp);
      td := t*d;
      Append(res,res_last,td);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end "*";

  procedure Apply ( L : in out List; t : in Transfo ) is

    res : constant List := t*L;

  begin
    Copy(res,L);
  end Apply;

  function Reduce ( L : List; i : integer32 ) return List is

    tmp,res,res_last : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      declare
        d1 : constant Link_to_Vector := Head_Of(tmp);
        d2 : constant Link_to_Vector := Reduce(d1,i);
      begin
       -- Append_Diff(res,res_last,d2);      -- be aware of duplicates !
        Append(res,res_last,d2);      -- be aware of duplicates !
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Reduce;

  procedure Reduce ( L : in out List; i : in integer32 ) is

    res : constant List := Reduce(L,i); 

  begin
    Copy(res,L);
  end Reduce;

  function Insert ( L : List; i,a : integer32 ) return List is

    tmp,res,res_last : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      declare
        d1 : constant Link_to_Vector := Head_Of(tmp);
        d2 : constant Link_to_Vector := Insert(d1,i,a);
      begin
        Append(res,res_last,d2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert;

  procedure Insert ( L : in out List; i,a : in integer32 ) is

    res : constant List := Insert(L,i,a);

  begin
    Deep_Clear(L);
    L := res;
  end Insert;

  function Transform_and_Reduce
             ( t : Transfo; i : integer32; L : List ) return List is

    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      declare
        d  : constant Link_to_Vector := Head_Of(tmp);
        td : constant Vector(d'range) := t*d.all;
        dr : constant Link_to_Vector := new Vector'(Reduce(td,i));
      begin
        Append(res,res_last,dr);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Transform_and_Reduce;

  procedure Transform_and_Reduce 
              ( t : in Transfo; i : in integer32; L : in out List ) is

    res : constant List := Transform_and_Reduce(t,i,L);

  begin
    Deep_Clear(L);
    L := res;
  end Transform_and_Reduce;

  function Insert_and_Transform
             ( L : List; i,a : integer32; t : Transfo ) return List is

    tmp,res,res_last : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      declare
        d : constant Link_to_Vector
          := Insert_and_Transform(Head_Of(tmp),i,a,t);
      begin
        Append(res,res_last,d);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert_and_Transform;

  procedure Insert_and_Transform
              ( L : in out List; i,a : in integer32; t : in Transfo ) is

    res : constant List := Insert_and_Transform(L,i,a,t);

  begin
    Deep_Clear(L);
    L := res;
  end Insert_and_Transform;

end Transforming_Integer32_Vector_Lists;
