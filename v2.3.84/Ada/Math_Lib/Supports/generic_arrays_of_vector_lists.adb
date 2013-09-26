with unchecked_deallocation;

package body Generic_Arrays_of_Vector_Lists is

-- CREATORS :

  function Deep_Create ( v : Array_of_VecVecs ) return Array_of_Lists is

    res : Array_of_Lists(v'range);

  begin
    for i in v'range loop
      if v(i) /= null
       then res(i) := Deep_Create(v(i).all);
      end if;   
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( v : Array_of_VecVecs ) return Array_of_Lists is

    res : Array_of_Lists(v'range);

  begin
    for i in v'range loop
      if v(i) /= null
       then res(i) := Shallow_Create(v(i).all);
      end if;   
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( L : Array_of_Lists ) return Array_of_VecVecs is

    res : Array_of_VecVecs(L'range);

  begin
    for i in res'range loop
      res(i) := new VecVec'(Deep_Create(L(i)));
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( L : Array_of_Lists ) return Array_of_VecVecs is

    res : Array_of_VecVecs(L'range);

  begin
    for i in res'range loop
      res(i) := new VecVec'(Shallow_Create(L(i)));
    end loop;
    return res;
  end Shallow_Create;

-- COMPARISON and COPYING :

  function Equal ( L1,L2 : Array_of_Lists ) return boolean is
  begin
    if L1'first /= L2'first or else L1'last /= L2'last then
      return false;
    else
      for k in L1'range loop
        if not Equal(L1(k),L2(k))
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Equal;

  function Equal ( L1,L2 : Link_to_Array_of_Lists ) return boolean is
  begin
    if L1 = null and then L2 /= null then
      return false;
    elsif L2 = null then
      return true;
    else
      return Equal(L1.all,L2.all);
    end if;
  end Equal;

  procedure Copy ( L1 : in Array_of_Lists; L2 : in out Array_of_Lists ) is
  begin
    for k in L1'range loop
      Copy(L1(k),L2(k));
    end loop;
  end Copy;

-- SELECTOR :

  function Length_Of ( L : Array_of_Lists ) return natural32 is

    res : natural32 := 0;
 
  begin
    for i in L'range loop
      res := res + Length_Of(L(i));
    end loop;
    return res;
  end Length_Of;

-- DESTRUCTORS :

  procedure free is
    new unchecked_deallocation(Array_of_Lists,Link_to_Array_of_Lists);

  procedure Deep_Clear ( L : in out Array_of_Lists ) is
  begin
    for k in L'range loop
      Deep_Clear(L(k));
    end loop;
  end Deep_Clear;

  procedure Shallow_Clear ( L : in out Array_of_Lists ) is
  begin
    for k in L'range loop
      Shallow_Clear(L(k));
    end loop;
  end Shallow_Clear;

  procedure Deep_Clear ( L : in out Link_to_Array_of_Lists ) is
  begin
    if L /= null
     then Deep_Clear(L.all); free(L);
    end if;
  end Deep_Clear;

  procedure Shallow_Clear ( L : in out Link_to_Array_of_Lists ) is
  begin
    if L /= null
     then Shallow_Clear(L.all); free(L);
    end if;
  end Shallow_Clear;

end Generic_Arrays_of_Vector_Lists;
