with unchecked_deallocation;

package body Generic_Polynomial_Vectors is

-- EVALUATORS :

  function Eval ( p : Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector is

    res : Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomials.Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Link_to_Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector is

    res : Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomials.Eval(p(i),x);
    end loop;
    return res;
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( p : in out Polynomial_Vector ) is
  begin
    for i in p'range loop
      Polynomials.Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Polynomial_Vector ) is

    procedure free is
      new unchecked_deallocation(Polynomial_Vector,Link_to_Polynomial_Vector);

  begin
    if p /= null then
      Clear(p.all);
      free(p);
    end if;
  end Clear;

end Generic_Polynomial_Vectors;
