package body Integer32_Vectors_Utilities is

  function Pivot ( v : Vector ) return integer32 is
  begin
    for i in v'range loop
      if v(i) /= 0
       then return i;
      end if;
    end loop;
    return (v'last + 1);
  end Pivot;

  function Pivot ( v : Link_to_Vector ) return integer32 is
  begin
    if v = null
     then return 0;
     else return Pivot(v.all);
    end if;
  end Pivot;

  function Reduce ( v : Vector; i : integer32 ) return Vector is

    res : Vector(v'first..(v'last-1));

  begin
    for j in res'range loop
      if j < i
       then res(j) := v(j);
       else res(j) := v(j+1);
      end if;
    end loop;
    return res;
  end Reduce;

  function Reduce ( v : Link_to_Vector; i : integer32 )
                  return Link_to_Vector is
  begin
    if v = null then
      return v;
    else
      declare
        res : constant Link_to_Vector := new Vector'(Reduce(v.all,i));
      begin
        return res;
      end;
    end if;
  end Reduce;

  procedure Reduce ( v : in out Link_to_Vector; i : in integer32 ) is
  begin
    if v /= null then
      declare
        res : constant Vector := Reduce(v.all,i);
      begin
        Clear(v);
        v := new Vector'(res);
      end;
    end if; 
  end Reduce;

  function Insert ( v : Vector; i,a : integer32 ) return Vector is

    res : Vector(v'first..(v'last+1));

  begin
    for j in res'first..(i-1) loop
      res(j) := v(j);
    end loop;
    res(i) := a;
    for j in (i+1)..res'last loop
      res(j) := v(j-1);
    end loop;
    return res;
  end Insert;

  function Insert ( v : Link_to_Vector; i,a : integer32 )
                  return Link_to_Vector is

    res : Link_to_Vector;

  begin
    if v = null
     then res := new Vector'(i..i => a);
     else res := new Vector'(Insert(v.all,i,a));
    end if;
    return res;
  end Insert;

  procedure Insert ( v : in out Link_to_Vector; i,a : in integer32 ) is
  begin
    if v /= null then
      declare
        res : constant Vector := Insert(v.all,i,a);
      begin
        Clear(v);
        v := new Vector'(res);
      end;
    end if;
  end Insert;

  function Insert_and_Transform
             ( v : Vector; i,a : integer32; t : Transfo ) return Vector is

    res : Vector(v'first..v'last+1) := Insert(v,i,a);

  begin
    Apply(t,res);
    return res;
  end Insert_and_Transform;

  procedure Insert_and_Transform
              ( v : in out Link_to_Vector; i,a : in integer32;
                t : in Transfo ) is

    res : Link_to_Vector;

  begin
    res := Insert_and_Transform(v,i,a,t);
    Clear(v);
    v := res;
  end Insert_and_Transform;

  function Insert_and_Transform
             ( v : Link_to_Vector; i,a : integer32; t : Transfo )
             return Link_to_Vector is

    res : Link_to_Vector;

  begin
    if v = null
     then res := Insert(v,i,a); Apply(t,res.all);
     else res := new Vector'(Insert_and_Transform(v.all,i,a,t));
    end if;
    return res;
  end Insert_and_Transform; 

end Integer32_Vectors_Utilities;
