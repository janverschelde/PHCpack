with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package body Parallel_Directions is

-- AUXILIARY FUNCTIONS :

  function Pivot ( v : Standard_Integer_Vectors.Vector ) return integer is
  begin
    for i in reverse v'range loop
      if v(i) /= 0
       then return i;
      end if;
    end loop;
    return v'first-1;
  end Pivot;

  function Normalize ( v : Standard_Integer_Vectors.Vector; p : integer )
                     return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(v'range);
    dvp : constant double_float := double_float(v(p));

  begin
    for i in res'first..p-1 loop
      res(i) := double_float(v(i))/dvp;
    end loop;
    res(p) := 1.0;
    for i in p+1..res'last loop
      res(i) := 0.0;
    end loop;
    return res;
  end Normalize;

  function ">" ( v,w : Standard_Floating_Vectors.Link_to_Vector )
               return boolean is
  begin
    for i in v'range loop
      if abs(v(i) - w(i)) > 1.0E-8 then
        if v(i) > w(i)
         then return true;
         else return false;
        end if;
      end if;
    end loop;
    return false;
  end ">";

  function Equal ( v,w : Standard_Floating_Vectors.Link_to_Vector )
                 return boolean is
  begin
    for i in v'range loop
      if abs(v(i) - w(i)) > 1.0E-8
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  procedure Swap ( v : in out Standard_Natural_Vectors.Vector;
                   i,j : in integer ) is

    temp : constant natural := v(i);

  begin
    v(i) := v(j);
    v(j) := temp;
  end Swap;

  function Min ( v : Standard_Floating_VecVecs.Array_of_VecVecs;
                 sv : Standard_Natural_VecVecs.VecVec;
                 ind : Standard_Natural_Vectors.Vector ) return integer is

    res,k : integer;
    p : integer := v'first;
    min : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while ind(p) > v(p)'last loop
      p := p + 1;
      exit when p > ind'last;
    end loop;
    if p > ind'last then
      return -1;
    else
      k := sv(p)(ind(p));
      min := v(p)(k);
      res := p;
      for i in p+1..v'last loop
        if ind(i) <= v(i)'last then
          k := sv(i)(ind(i));
          if min > v(i)(k) then
            res := i;
            min := v(i)(k);
          end if;
        end if;
      end loop;
      return res;
    end if;
  end Min;

  function Last_Index ( v : Standard_Floating_VecVecs.Array_of_VecVecs;
                        ind : Standard_Natural_Vectors.Vector )
                      return integer is

    res : integer := -1;

  begin
    for i in v'range loop
      if ind(i) <= v(i)'last then
        if res = -1
         then res := i;
         else return -1;
        end if;
      end if;
    end loop;
    if res = -1
     then return 0;
     else return res;
    end if;
  end Last_Index;

end Parallel_Directions;
