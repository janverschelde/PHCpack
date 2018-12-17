with unchecked_deallocation;

package body Generic_Monomial_Vectors is

-- EVALUATORS :

  function Eval ( v : Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := Monomials.Eval(v(v'first),x);

    use Ring;

  begin
    for i in v'first+1..v'last loop
      res := res + Monomials.Eval(v(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( v : Link_to_Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := Monomials.Eval(v(v'first),x);

    use Ring;

  begin
    for i in v'first+1..v'last loop
      res := res + Monomials.Eval(v(i),x);
    end loop;
    return res;
  end Eval;

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

end Generic_Monomial_Vectors;
