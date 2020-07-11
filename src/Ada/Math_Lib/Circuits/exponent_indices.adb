package body Exponent_Indices is

  function Index_Size ( v : Standard_Integer_Vectors.Vector ) 
                      return integer32 is

    res : integer32 := 0;

  begin
    for i in v'range loop
      if v(i) > 0
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Index_Size;

  function Index_Size ( v : Standard_Integer_Vectors.Link_to_Vector ) 
                      return integer32 is
  begin
    return Index_Size(v.all);
  end Index_Size;

  function Factor_Size ( v : Standard_Integer_Vectors.Link_to_Vector ) 
                       return integer32 is

    res : integer32 := 0;

  begin
    for i in v'range loop
      if v(i) > 1
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Factor_Size;

  function Exponent_Index
             ( xp : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..Index_Size(xp));
    idx : integer32 := 0;

  begin
    for k in xp'range loop
      if xp(k) > 0 then
        idx := idx + 1;
        res(idx) := k;
      end if;
    end loop;
    return res;
  end Exponent_Index;

  function Exponent_Index
             ( xp : Standard_Integer_Vectors.Link_to_Vector )
             return Standard_Integer_Vectors.Link_to_Vector is

    indices : constant Standard_Integer_Vectors.Vector
            := Exponent_Index(xp.all);
    res : constant Standard_Integer_Vectors.Link_to_Vector
        := new Standard_Integer_Vectors.Vector'(indices);

  begin
    return res;
  end Exponent_Index;

  function Factor_Index
             ( xp : Standard_Integer_Vectors.Link_to_Vector )
             return Standard_Integer_Vectors.Link_to_Vector is

    res : Standard_Integer_Vectors.Link_to_Vector := null;
    fsz : constant integer32 := Factor_Size(xp);

  begin
    if fsz > 0 then
      declare
        fac : Standard_Integer_Vectors.Vector(1..fsz);
        idx : integer32 := 0;
      begin
        for k in xp'range loop
          if xp(k) > 1 then
            idx := idx + 1;
            fac(idx) := k;
          end if;
        end loop;
        res := new Standard_Integer_Vectors.Vector'(fac);
      end;
    end if;
    return res;
  end Factor_Index;

  function Exponent_Index
             ( xp : Standard_Integer_VecVecs.VecVec )
             return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(xp'range);

  begin
    for i in xp'range loop
      res(i) := Exponent_Index(xp(i));
    end loop;
    return res;
  end Exponent_Index;

  function Factor_Index
             ( xp : Standard_Integer_VecVecs.VecVec )
             return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(xp'range);

  begin
    for i in xp'range loop
      res(i) := Factor_Index(xp(i));
    end loop;
    return res;
  end Factor_Index;

  function Maxima ( xp : Standard_Integer_VecVecs.VecVec )
                  return Standard_Integer_Vectors.Vector is

    lpt : Standard_Integer_Vectors.Link_to_Vector := xp(xp'first);
    res : Standard_Integer_Vectors.Vector(lpt'range) := lpt.all;

  begin
    for k in xp'first+1..xp'last loop
      lpt := xp(k);
      for i in res'range loop
        if lpt(i) > res(i)
         then res(i) := lpt(i);
        end if;
      end loop;
    end loop;
    return res;
  end Maxima;

  function Polynomial_Degree
             ( xp : Standard_Integer_Vectors.Link_to_Vector )
             return integer32 is

    res : integer32 := 0;

    use Standard_Integer_Vectors;

  begin
    if xp = null then
      return -1;
    else
      for k in xp'range loop
        res := res + xp(k);
      end loop;
      return res;
    end if;
  end Polynomial_Degree;

  function Polynomial_Degree
             ( xp : Standard_Integer_VecVecs.VecVec ) return integer32 is

    res : integer32 := -1;
    deg : integer32;

  begin
    for k in xp'range loop
      deg := Polynomial_Degree(xp(k));
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Polynomial_Degree;

end Exponent_Indices;
