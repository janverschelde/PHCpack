with Standard_Complex_Numbers;
with Standard_Mathematical_Functions;
with DoblDobl_Complex_Numbers;
with DoblDobl_Mathematical_Functions;
with QuadDobl_Complex_Numbers;
with QuadDobl_Mathematical_Functions;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
--with Floating_Lifting_Utilities;         use Floating_Lifting_Utilities;
--with Floating_Integer_Convertors;        use Floating_Integer_Convertors;

package body Polyhedral_Coefficient_Homotopies is

  outside_upper : constant integer32 := 1579;

  function Minimum ( v : Standard_Integer_VecVecs.VecVec ) return integer32 is

    min : integer32 := 0;
    tmp : integer32;

  begin
    for i in v'range loop
      declare
        vv : constant Standard_Integer_Vectors.Vector := v(i).all;
      begin
        for j in vv'range loop
          if vv(j) /= 0 then 
            if vv(j) < 0
             then tmp := -vv(j);
             else tmp := vv(j);
            end if;
            if (min = 0) or else (tmp < min) then
              min := tmp;
            end if;
          end if;
        end loop;
      end;
    end loop;
    return min;
  end Minimum;

  function Minimum ( v : Standard_Floating_VecVecs.VecVec )
                   return double_float is

    tol : constant double_float := 1.0E-8;
    min : double_float := 0.0;
    tmp : double_float;

  begin
    for i in v'range loop
      declare
        vv : constant Standard_Floating_Vectors.Vector := v(i).all;
      begin
        for j in vv'range loop
          tmp := abs(vv(j));
          if tmp > tol then
            if (min = 0.0) or else (tmp < min)
             then min := tmp;
            end if;
          end if;
        end loop;
      end;
    end loop;
    return min;
  end Minimum;

  function Scale ( v : Standard_Integer_Vectors.Vector; s : integer32 )
                 return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(v'range);

  begin
    if (s = 0) or (s = 1) then
      for i in res'range loop
        res(i) := double_float(v(i));
      end loop;
    else
      for i in res'range loop
        res(i) := double_float(v(i))/double_float(s);
      end loop;
    end if;
    return res;
  end Scale;

  function Scale ( v : Standard_Floating_Vectors.Vector; s : double_float )
                 return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(v'range);
    tol : constant double_float := 1.0E-8;

  begin
    if (abs(s) < tol) or (s = 1.0) then
      for i in res'range loop
        res(i) := v(i);
      end loop;
    else
      for i in res'range loop
        res(i) := v(i)/s;
      end loop;
    end if;
    return res;
  end Scale;

  procedure Scale ( v : in out Standard_Floating_Vectors.Vector;
                    s : in double_float ) is

    tol : constant double_float := 1.0E-8;

  begin
    if (abs(s) > tol) and (s /= 1.0) then
      for i in v'range loop
        v(i) := v(i)/s;
      end loop;
    end if;
  end Scale;

  function Scale ( v : Standard_Integer_VecVecs.VecVec )
                 return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(v'range);
    min : constant integer32 := Minimum(v);

  begin
    for i in v'range loop
      declare
        sv : constant Standard_Floating_Vectors.Vector
           := Scale(v(i).all,min);
      begin
        res(i) := new Standard_Floating_Vectors.Vector(sv'range);
        for j in sv'range loop
          res(i)(j) := sv(j);
        end loop;
      end;  -- this detour was set up for GNAT 3.07 ...
     -- res(i) := new Standard_Floating_Vectors.Vector'(Scale(v(i).all));
    end loop;
    return res;
  end Scale;

  function Scale ( v : Standard_Floating_VecVecs.VecVec )
                 return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(v'range);
    min : constant double_float := Minimum(v);

  begin
    for i in v'range loop
      declare
        sv : constant Standard_Floating_Vectors.Vector
           := Scale(v(i).all,min);
      begin
        res(i) := new Standard_Floating_Vectors.Vector(sv'range);
        for j in sv'range loop
          res(i)(j) := sv(j);
        end loop;
      end;  -- this detour was set up for GNAT 3.07 ...
     -- res(i) := new Standard_Floating_Vectors.Vector'(Scale(v(i).all));
    end loop;
    return res;
  end Scale;

  procedure Scale ( v : in out Standard_Floating_VecVecs.VecVec ) is

    min : constant double_float := Minimum(v);

  begin
    for i in v'range loop
      Scale(v(i).all,min);
    end loop;
  end Scale;

  procedure Shift ( v : in out Standard_Integer_Vectors.Vector ) is

    min : integer32 := v(v'first);

  begin
    for i in v'first+1..v'last loop
      if v(i) < min
       then min := v(i);
      end if;
    end loop;
    if min /= 0 then
      for i in v'range loop
        v(i) := v(i) - min;
      end loop;
    end if;
  end Shift;

  procedure Shift ( v : in out Standard_Floating_Vectors.Vector ) is

    min : double_float := v(v'first);

  begin
    for i in v'first+1..v'last loop
      if v(i) < min
       then min := v(i);
      end if;
    end loop;
    if min /= 0.0 then
      for i in v'range loop
        v(i) := v(i) - min;
      end loop;
    end if;
  end Shift;

  procedure Search_Lifting
              ( L : in Lists_of_Floating_Vectors.List;
                x : in Standard_Integer_Vectors.Vector;
                found : out boolean; y : out double_float ) is

    tmp : Lists_of_Floating_Vectors.List := L;
    lpt : Standard_Floating_Vectors.Link_to_Vector;
 
  begin
    found := false;
    while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
      lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
      found := true;
      for i in x'range loop
        if integer32(lpt(i)) /= x(i)
         then found := false;
        end if;
        exit when not found;
      end loop;
      if found 
       then y := lpt(lpt'last);
      end if;
      exit when found;
      tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
    end loop;
  end Search_Lifting;

  function Power_Transform
              ( e : Standard_Integer_VecVecs.VecVec;
                s : Lists_of_Integer_Vectors.List;
                normal : Standard_Integer_Vectors.Vector ) 
              return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(e'range);
    lei : Standard_Integer_Vectors.Vector(normal'first..normal'last-1);
    found : boolean;
    lif : integer32;
    upper : constant integer32 := outside_upper; -- for outside points

    use Standard_Integer_Vectors;

  begin
    for i in e'range loop
      lei := e(i).all;
      Search_Lifting(s,lei,found,lif);
      if not found
       then res(i) := upper;
       else res(i) := lei*normal(lei'range) + normal(normal'last)*lif;
      end if;
    end loop;
    Shift(res);
    return res;
  end Power_Transform;

  function Power_Transform
              ( e : Standard_Integer_VecVecs.VecVec;
                s : Lists_of_Floating_Vectors.List;
                normal : Standard_Floating_Vectors.Vector ) 
              return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(e'range);
    found : boolean;
    lif : double_float;
    upper : constant double_float
          := double_float(outside_upper);  -- for outside points

  begin
    for i in e'range loop
      Search_Lifting(s,e(i).all,found,lif);
      if not found then
        res(i) := upper;
      else
       -- res(i) := lei*normal(lei'range) + normal(normal'last)*lif;
        res(i) := normal(normal'last)*lif;
        for j in e(i)'range loop
          res(i) := res(i) + normal(j)*double_float(e(i)(j));
        end loop;
      end if;
    end loop;
    Shift(res);
    return res;
  end Power_Transform;

  procedure Power_Transform
              ( e : in Standard_Integer_VecVecs.VecVec;
                s : in Lists_of_Floating_Vectors.List;
                normal : in Standard_Floating_Vectors.Vector;
                L : out Standard_Floating_Vectors.Vector ) is

    found : boolean;
    lif : double_float;
    upper : constant double_float
          := double_float(outside_upper);  -- for outside points

  begin
    for i in e'range loop
      Search_Lifting(s,e(i).all,found,lif);
      if not found then
        L(i) := upper;
      else
        L(i) := normal(normal'last)*lif;
        for j in e(i)'range loop
          L(i) := L(i) + normal(j)*double_float(e(i)(j));
        end loop;
      end if;
    end loop;
    Shift(L);
  end Power_Transform;

  function Power_Transform
              ( e : Exponent_Vectors.Exponent_Vectors_Array;
                s : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix,normal : Standard_Integer_Vectors.Vector )
              return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(e'range);
    cnt : integer32 := res'first;
 
  begin
    for i in mix'range loop
      declare
        rescnt : constant Standard_Integer_Vectors.Vector
               := Power_Transform(e(cnt).all,s(i),normal);
      begin
        res(cnt) := new Standard_Integer_Vectors.Vector(rescnt'range);
        for j in rescnt'range loop
          res(cnt)(j) := rescnt(j);
        end loop;      
      end;
     -- Shift(res(cnt).all);
      for k in 1..(mix(i)-1) loop
        res(cnt+k) := new Standard_Integer_Vectors.Vector(res(cnt)'range);
        for j in res(cnt)'range loop
          res(cnt+k)(j) := res(cnt)(j);
        end loop;
      end loop;
      cnt := cnt + mix(i);
    end loop;
    return res;
  end Power_Transform;

  function Power_Transform
              ( e : Exponent_Vectors.Exponent_Vectors_Array;
                s : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mix : Standard_Integer_Vectors.Vector;
                normal : Standard_Floating_Vectors.Vector )
              return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(e'range);
    cnt : integer32 := res'first;
 
  begin
    for i in mix'range loop
      declare
        rescnt : constant Standard_Floating_Vectors.Vector
               := Power_Transform(e(cnt).all,s(i),normal);
      begin
        res(cnt) := new Standard_Floating_Vectors.Vector(rescnt'range);
        for j in rescnt'range loop
          res(cnt)(j) := rescnt(j);
        end loop;      
      end;
     -- Shift(res(cnt).all);
      for k in 1..(mix(i)-1) loop
        res(cnt+k) := new Standard_Floating_Vectors.Vector(res(cnt)'range);
        for j in res(cnt)'range loop
          res(cnt+k)(j) := res(cnt)(j);
        end loop;
      end loop;
      cnt := cnt + mix(i);
    end loop;
    return res;
  end Power_Transform;

  procedure Power_Transform
              ( e : in Exponent_Vectors.Exponent_Vectors_Array;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Vector;
                normal : in Standard_Floating_Vectors.Vector;
                L : in out Standard_Floating_VecVecs.VecVec ) is

    cnt : integer32 := L'first;
 
  begin
    for i in mix'range loop
      Power_Transform(e(cnt).all,s(i),normal,L(cnt).all);
     -- Shift(L(cnt).all);
      for k in 1..(mix(i)-1) loop
        for j in L(cnt)'range loop
          L(cnt+k)(j) := L(cnt)(j);
        end loop;
      end loop;
      cnt := cnt + mix(i);
    end loop;
    Scale(L);
  end Power_Transform;

  procedure Eval ( c : in Standard_Complex_Vectors.Vector;
                   t : in double_float;
                   m : in Standard_Integer_Vectors.Vector;
                   ctm : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

  begin
    for i in ctm'range loop
      ctm(i) := c(i)*Create((t**integer(m(i))));
    end loop;
  end Eval;

  procedure Eval ( c : in Standard_Complex_Vectors.Vector;
                   t : in double_float;
                   m : in Standard_Floating_Vectors.Vector;
                   ctm : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers,Standard_Mathematical_Functions;

  begin
    for i in ctm'range loop
      if (REAL_PART(c(i)) = 0.0) and (IMAG_PART(c(i)) = 0.0)
       then ctm(i) := Create(0.0);
       else ctm(i) := c(i)*Create(t**m(i));
      end if;
    end loop;
  end Eval;

  procedure Eval ( c : in Standard_Complex_VecVecs.VecVec;
                   t : in double_float; 
                   m : in Standard_Integer_VecVecs.VecVec;
                   ctm : in out Standard_Complex_VecVecs.VecVec ) is
  begin
    for i in ctm'range loop
      Eval(c(i).all,t,m(i).all,ctm(i).all);
    end loop;
  end Eval;

  procedure Eval ( c : in Standard_Complex_VecVecs.VecVec;
                   t : in double_float;
                   m : in Standard_Floating_VecVecs.VecVec;
                   ctm : in out Standard_Complex_VecVecs.VecVec ) is
  begin
    for i in ctm'range loop
      Eval(c(i).all,t,m(i).all,ctm(i).all);
    end loop;
  end Eval;

  procedure Eval ( c : in DoblDobl_Complex_Vectors.Vector;
                   t : in double_double;
                   m : in Standard_Integer_Vectors.Vector;
                   ctm : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

  begin
    for i in ctm'range loop
      ctm(i) := c(i)*Create((t**integer(m(i))));
    end loop;
  end Eval;

  --function PowerFloat ( t : double_double;
  --                      m : double_float ) return double_double is

  -- DESCRIPTION :
  --   The approximation hi_part(t)**m is better than t**integer(m).

  --  use Standard_Mathematical_Functions;

  --  tsd : double_float := hi_part(t);
  --  rsd : double_float := tsd**m;
  --  res : double_double := create(rsd);

  --begin
  --  return res;
  --end PowerFloat;

  --function PowerFloat ( t : quad_double;
  --                      m : double_float ) return quad_double is

  -- DESCRIPTION :
  --   The approximation hi_part(t)**m is better than t**integer(m).

  --  use Standard_Mathematical_Functions;

  --  tsd : double_float := hihi_part(t);
  --  rsd : double_float := tsd**m;
  --  res : quad_double := create(rsd);

  --begin
  --  return res;
  --end PowerFloat;

  procedure Eval ( c : in DoblDobl_Complex_Vectors.Vector;
                   t : in double_double;
                   m : in Standard_Floating_Vectors.Vector;
                   ctm : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers,DoblDobl_Mathematical_Functions;

    zero : constant double_double := create(0.0);
    tmi : double_double;

  begin
    for i in ctm'range loop
      if (REAL_PART(c(i)) = zero) and (IMAG_PART(c(i)) = zero) then
        ctm(i) := Create(zero);
      else
        tmi := t**m(i); -- PowerFloat(t,m(i)); -- t**integer(m(i));
        ctm(i) := c(i)*Create(tmi);
      end if;
    end loop;
  end Eval;

  procedure Eval ( c : in DoblDobl_Complex_VecVecs.VecVec;
                   t : in double_double; 
                   m : in Standard_Integer_VecVecs.VecVec;
                   ctm : in out DoblDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in ctm'range loop
      Eval(c(i).all,t,m(i).all,ctm(i).all);
    end loop;
  end Eval;

  procedure Eval ( c : in DoblDobl_Complex_VecVecs.VecVec;
                   t : in double_double;
                   m : in Standard_Floating_VecVecs.VecVec;
                   ctm : in out DoblDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in ctm'range loop
      Eval(c(i).all,t,m(i).all,ctm(i).all);
    end loop;
  end Eval;

  procedure Eval ( c : in QuadDobl_Complex_Vectors.Vector;
                   t : in quad_double;
                   m : in Standard_Integer_Vectors.Vector;
                   ctm : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

  begin
    for i in ctm'range loop
      ctm(i) := c(i)*Create((t**integer(m(i))));
    end loop;
  end Eval;

  procedure Eval ( c : in QuadDobl_Complex_Vectors.Vector;
                   t : in quad_double;
                   m : in Standard_Floating_Vectors.Vector;
                   ctm : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers,QuadDobl_Mathematical_Functions;
    zero : constant quad_double := create(0.0);
    tmi : quad_double;

  begin
    for i in ctm'range loop
      if (REAL_PART(c(i)) = zero) and (IMAG_PART(c(i)) = zero) then
        ctm(i) := Create(zero);
      else
        tmi := t**m(i); -- PowerFloat(t,m(i));
        ctm(i) := c(i)*Create(tmi);
      end if;
    end loop;
  end Eval;

  procedure Eval ( c : in QuadDobl_Complex_VecVecs.VecVec;
                   t : in quad_double; 
                   m : in Standard_Integer_VecVecs.VecVec;
                   ctm : in out QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in ctm'range loop
      Eval(c(i).all,t,m(i).all,ctm(i).all);
    end loop;
  end Eval;

  procedure Eval ( c : in QuadDobl_Complex_VecVecs.VecVec;
                   t : in quad_double;
                   m : in Standard_Floating_VecVecs.VecVec;
                   ctm : in out QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in ctm'range loop
      Eval(c(i).all,t,m(i).all,ctm(i).all);
    end loop;
  end Eval;

end Polyhedral_Coefficient_Homotopies;
