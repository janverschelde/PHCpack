with unchecked_deallocation;

package body Generic_Monomial_Vectors is

-- SELECTORS :

  function Degree ( v : Monomial_Vector ) return integer32 is

    res : integer32 := Monomials.Degree(v(v'first));
    deg : integer32;

  begin
    for i in v'first+1..v'last loop
      deg := Monomials.Degree(v(i));
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Degree;

  function Degree ( v : Link_to_Monomial_Vector ) return integer32 is
  begin
    if v = null
     then return -1;
     else return Degree(v.all);
    end if;
  end Degree;

  function Degree ( p : Polynomial ) return integer32 is

    res : integer32 := Degree(p.mons);

  begin
    if res >= 0 then
      return res;
    elsif not Ring.Equal(p.cff0,Ring.zero) then
      return 0;
    else
      return -1;
    end if;
  end Degree;

  function Degree ( p : Link_to_Polynomial ) return integer32 is
  begin
    if p = null
     then return -1;
     else return Degree(p.all);
    end if;
  end Degree;

-- EVALUATION and DIFFERENTIATION :

  function Eval ( v : Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := Monomials.Eval(v(v'first),x);
    eva : Ring.number;

    use Ring;

  begin
    for i in v'first+1..v'last loop
      eva := Monomials.Eval(v(i),x);
      Add(res,eva);
      Clear(eva);
    end loop;
    return res;
  end Eval;

  function Eval ( v : Link_to_Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := Monomials.Eval(v(v'first),x);
    eva : Ring.number;

    use Ring;

  begin
    for i in v'first+1..v'last loop
      eva := Monomials.Eval(v(i),x);
      Add(res,eva);
      Clear(eva);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Polynomial;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := Eval(p.mons,x);

  begin
    Ring.add(res,p.cff0);
    return res;
  end Eval;

  function Eval ( p : Link_to_Polynomial;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number;

  begin
    if p = null
     then Ring.copy(Ring.zero,res);
     else res := Eval(p.all,x);
    end if;
    return res;
  end Eval;

  procedure Diff ( v : in Monomial_Vector; x : in Vectors.Vector;
                   yd,wrk : in out Vectors.Vector ) is
  begin
    for i in wrk'range loop
      Ring.copy(Ring.zero,yd(i));
    end loop;
    for i in v'range loop
      Monomials.Diff(v(i),x,wrk);
      for j in 1..integer32(v(i).nvr) loop
        Ring.add(yd(integer32(v(i).pos(j))),wrk(j));
      end loop;
    end loop;
  end Diff;

  procedure Diff ( v : in Monomial_Vector; x : in Vectors.Vector;
                   yd : in out Vectors.Vector ) is

    wrk : Vectors.Vector(x'range);

  begin
    Diff(v,x,yd,wrk);
  end Diff;

  procedure Diff ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                   yd,wrk : in out Vectors.Vector ) is
  begin
    if v /= null
     then Diff(v.all,x,yd,wrk);
    end if;
  end Diff;

  procedure Diff ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                   yd : in out Vectors.Vector ) is
  begin
    if v /= null
     then Diff(v.all,x,yd);
    end if;
  end Diff;

  procedure Speel ( v : in Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is

    mon : Monomials.Link_to_Monomial;
    val : Ring.number;

  begin
    for i in wrk'range loop
      Ring.copy(Ring.zero,yd(i));
    end loop;
    mon := v(v'first);
    Monomials.Speel(mon,x,y,wrk);
    for j in 1..integer32(mon.nvr) loop
      Ring.add(yd(integer32(mon.pos(j))),wrk(j));
    end loop;
    for i in v'first+1..v'last loop
      mon := v(i);
      Monomials.Speel(mon,x,val,wrk);
      Ring.add(y,val);
      for j in 1..integer32(mon.nvr) loop
        Ring.add(yd(integer32(mon.pos(j))),wrk(j));
      end loop;
    end loop;
  end Speel;

  procedure Speel ( v : in Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    wrk : Vectors.Vector(x'range);

  begin
    Speel(v,x,y,yd,wrk);
  end Speel;

  procedure Speel ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    if v /= null
     then Speel(v.all,x,y,yd,wrk);
    end if;
  end Speel;

  procedure Speel ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is
  begin
    if v /= null
     then Speel(v.all,x,y,yd);
    end if;
  end Speel;

  procedure Speel ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    Speel(p.mons,x,y,yd,wrk);
    Ring.add(y,p.cff0);
  end Speel;

  procedure Speel ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is
  begin
    Speel(p.mons,x,y,yd);
    Ring.add(y,p.cff0);
  end Speel;

  procedure Speel ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    if p /= null
     then Speel(p.all,x,y,yd,wrk);
    end if;
  end Speel;

  procedure Speel ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is
  begin
    if p /= null
     then Speel(p.all,x,y,yd);
    end if;
  end Speel;

-- DESTRUCTORS :

  procedure Clear ( v : in out Monomial_Vector ) is
  begin
    for i in v'range loop
      Monomials.Clear(v(i));
    end loop;
  end Clear;

  procedure Clear ( v : in out Link_to_Monomial_Vector ) is

    procedure free is
      new unchecked_deallocation(Monomial_Vector,Link_to_Monomial_Vector);

  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

  procedure Clear ( p : in out Polynomial ) is
  begin
    Clear(p.mons);
  end Clear;

  procedure Clear ( p : in out Link_to_Polynomial ) is

    procedure free is
      new unchecked_deallocation(Polynomial,Link_to_Polynomial);

  begin
    if p /= null then
      Clear(p.all);
      free(p);
    end if;
  end Clear;

end Generic_Monomial_Vectors;
