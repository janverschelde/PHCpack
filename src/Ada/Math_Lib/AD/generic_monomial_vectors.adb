with unchecked_deallocation;

package body Generic_Monomial_Vectors is

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

  begin
    if p = null
     then return Ring.zero;
     else return Eval(p.all,x);
    end if;
  end Eval;

  procedure Speel ( v : in Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    for i in wrk'range loop
      wrk(i) := Ring.zero;
    end loop;
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
